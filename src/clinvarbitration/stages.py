from datetime import datetime
from functools import cache
from os.path import join
from typing import TYPE_CHECKING

from cpg_utils import Path, to_path
from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import get_batch

from cpg_flow.stage import MultiCohortStage, stage

from clinvarbitration.jobs.AnnotateClinvarSnvsWithBcftools import annotate_clinvar_snvs
from clinvarbitration.jobs.CopyLatestClinvarFiles import copy_latest_files
from clinvarbitration.jobs.GenerateNewClinvarSummary import generate_new_summary
from clinvarbitration.jobs.PackageForRelease import package_data_for_release
from clinvarbitration.jobs.Pm5TableGeneration import generate_pm5_data


if TYPE_CHECKING:
    from cpg_flow.stage import StageInput, StageOutput
    from cpg_utils import Path
    from cpg_flow.targets.multicohort import MultiCohort


@cache
def get_output_folder():
    """
    get the folder to use for this run
    sits in the cpg-common bucket, in a folder named YY-MM
    """

    return to_path(
        join(
            config_retrieve(['storage', 'common', 'analysis']),
            'clinvarbitration',
            datetime.now().strftime('%y-%m'),  # noqa: DTZ005
        ),
    )


def populate_job_meta(output_file: str):
    """
    populates some metadata for the job

    Args:
        output_file (str): path to the output file

    Returns:
        dict of metadata
    """

    from clinvarbitration import __version__ as clinvarbitration_version

    print(f'Generating output meta for {output_file}')
    return {'image': config_retrieve(['workflow', 'driver_image']), 'version': clinvarbitration_version}


@stage
class CopyLatestClinvarFiles(MultiCohortStage):
    def expected_outputs(self, mc: 'MultiCohort') -> 'dict[str, Path]':
        return {
            'submission_file': get_output_folder() / 'submission_summary.txt.gz',
            'variant_file': get_output_folder() / 'variant_summary.txt.gz',
        }

    def queue_jobs(self, mc: 'MultiCohort', inputs: 'StageInput') -> 'StageOutput':
        """
        run a wget copy of the relevant files into GCP
        """

        outputs = self.expected_outputs(mc)

        bash_job = copy_latest_files(submissions=outputs['submission_file'], variants=outputs['variant_file'])

        return self.make_outputs(data=outputs, jobs=bash_job, target=mc)


@stage(required_stages=[CopyLatestClinvarFiles], analysis_type='clinvarbitration', analysis_keys=['clinvar_decisions'])
class GenerateNewClinvarSummary(MultiCohortStage):
    def expected_outputs(self, mc: 'MultiCohort') -> 'dict[str, Path]':
        """
        a couple of files and a HT as Paths
        """
        return {
            'clinvar_decisions': get_output_folder() / 'clinvar_decisions.ht.tar.gz',
            'snv_vcf': get_output_folder() / 'clinvar_decisions.vcf.bgz',
        }

    def queue_jobs(self, mc: 'MultiCohort', inputs: 'StageInput') -> 'StageOutput':
        outputs = self.expected_outputs(mc)
        var_file = get_batch().read_input(inputs.as_str(mc, CopyLatestClinvarFiles, 'variant_file'))
        sub_file = get_batch().read_input(inputs.as_str(mc, CopyLatestClinvarFiles, 'submission_file'))
        job = generate_new_summary(
            var_file=var_file,
            sub_file=sub_file,
            output_root=str(outputs['clinvar_decisions']).removesuffix('.ht.tar.gz'),
        )
        return self.make_outputs(target=mc, data=outputs, jobs=job)


@stage(required_stages=[GenerateNewClinvarSummary])
class AnnotateClinvarSnvsWithBcftools(MultiCohortStage):
    """
    take the vcf output from the clinvar stage, and apply consequence annotations
    this uses BCFtools for speed and re-deployability, instead of VEP
    """

    def expected_outputs(self, mc: 'MultiCohort') -> 'Path':
        return get_output_folder() / 'clinvar_decisions.annotated.tsv'

    def queue_jobs(self, mc: 'MultiCohort', inputs: 'StageInput') -> 'StageOutput':
        output = self.expected_outputs(mc)
        snv_vcf = get_batch().read_input(inputs.as_str(mc, GenerateNewClinvarSummary, 'snv_vcf'))
        job = annotate_clinvar_snvs(snv_vcf=snv_vcf, output=output)
        return self.make_outputs(target=mc, jobs=job, data=output)


@stage(required_stages=[AnnotateClinvarSnvsWithBcftools])
class Pm5TableGeneration(MultiCohortStage):
    def expected_outputs(self, mc: 'MultiCohort') -> dict[str, 'Path']:
        """
        a single HT, compressed into a tarball
        """
        return {
            'pm5_ht': get_output_folder() / 'clinvar_decisions.pm5.ht.tar.gz',
            'pm5_json': get_output_folder() / 'clinvar_decisions.pm5.json',
        }

    def queue_jobs(self, mc: 'MultiCohort', inputs: 'StageInput') -> 'StageOutput':
        outputs = self.expected_outputs(mc)
        annotated_snvs = inputs.as_str(mc, AnnotateClinvarSnvsWithBcftools)

        job = generate_pm5_data(
            annotated_snvs=annotated_snvs,
            output_root=str(outputs['pm5_json']).removesuffix('.json'),
        )

        return self.make_outputs(target=mc, data=outputs, jobs=job)


@stage(
    required_stages=[
        GenerateNewClinvarSummary,
        AnnotateClinvarSnvsWithBcftools,
        Pm5TableGeneration,
    ],
    analysis_type='clinvarbitration',
    update_analysis_meta=populate_job_meta,
)
class PackageForRelease(MultiCohortStage):
    """
    takes the data created so far, and packages it up for release
    """

    def expected_outputs(self, multicohort: 'MultiCohort') -> 'Path':
        return get_output_folder() / 'clinvar_decisions.release.tar.gz'

    def queue_jobs(self, multicohort: 'MultiCohort', inputs: 'StageInput') -> 'StageOutput':
        """
        Localise all the previously generated data into a folder
        tarball it, and write out as a single file
        """
        output = self.expected_outputs(multicohort)
        pm5 = inputs.as_dict(multicohort, Pm5TableGeneration)
        clinvar_decisions = inputs.as_str(multicohort, GenerateNewClinvarSummary, 'clinvar_decisions')
        job = package_data_for_release(
            annotated_tsv=inputs.as_path(multicohort, AnnotateClinvarSnvsWithBcftools),
            pm5_json=pm5['pm5_json'],
            pm5_ht=pm5['pm5_ht'],
            clinvar_decisions=clinvar_decisions,
            output=output,
        )
        return self.make_outputs(multicohort, data=output, jobs=job)


@stage(
    analysis_keys=['release.tar.gz'],
    analysis_type='clinvarbitration',
    update_analysis_meta=populate_job_meta,
)
class ClinvarbitrationNextflow(MultiCohortStage):
    """
    Instead of us running this one way, and off-site users running it another way,
    this single stage executes the full NextFlow workflow
    We don't have to use this, but we should run it alongside our main workflow to ensure both are working
    """

    def expected_outputs(
        self,
        multicohort: 'MultiCohort',
    ) -> 'dict[str, Path]':
        return {
            'submission_raw.txt.gz': get_output_folder() / 'clinvar_decisions.submission_raw.txt.gz',
            'variant_raw.txt.gz': get_output_folder() / 'clinvar_decisions.variant_raw.txt.gz',
            'ht.tar.gz': get_output_folder() / 'clinvar_decisions.ht.tar.gz',
            'pm5.ht.tar.gz': get_output_folder() / 'clinvar_decisions.pm5.ht.tar.gz',
            'vcf.bgz': get_output_folder() / 'clinvar_decisions.vcf.bgz',
            'release.tar.gz': get_output_folder() / 'clinvar_decisions.release.tar.gz',
        }

    def queue_jobs(
        self,
        multicohort: 'MultiCohort',
        inputs: 'StageInput',
    ) -> 'StageOutput':
        """

        Args:
            multicohort ():
            inputs ():

        Returns:

        """

        outputs = self.expected_outputs(multicohort)

        # read in a reference genome
        ref_fa = get_batch().read_input(config_retrieve(['workflow', 'ref_fa']))

        # make a new job
        job = get_batch().new_job('Run ClinvArbitration Nextflow')

        job.image(config_retrieve(['workflow', 'driver_image']))

        # set some resource params
        job.storage('10Gi').memory('highmem').cpu(2)

        # set up one sprawling output group for all workflow results
        job.declare_resource_group(
            output={
                'submission_raw.txt.gz': '{root}/submission_summary.txt.gz',
                'variant_raw.txt.gz': '{root}/variant_summary.txt.gz',
                'clinvar_decisions.json': '{root}/clinvar_decisions.json',
                'ht.tar.gz': '{root}/clinvar_decisions.ht.tar.gz',
                'pm5.ht.tar.gz': '{root}/clinvar_decisions.pm5.ht.tar.gz',
                'pm5.json': '{root}/clinvar_decisions.pm5.json',
                'vcf.bgz': '{root}/clinvar_decisions.vcf.bgz',
                'vcf.bgz.tbi': '{root}/clinvar_decisions.vcf.bgz.tbi',
                'unfiltered.vcf.bgz': '{root}/clinvar_decisions.unfiltered.vcf.bgz',
                'unfiltered.vcf.bgz.tbi': '{root}/clinvar_decisions.unfiltered.vcf.bgz.tbi',
                'annotated.tsv': '{root}/clinvar_decisions.annotated.tsv',
                'release.tar.gz': '{root}/clinvar_decisions.release.tar.gz',
            },
        )

        gff3 = get_batch().read_input(config_retrieve(['references', 'ensembl_113', 'gff3']))

        # nextflow go brrrr
        job.command(
            f"""
            nextflow \
                -c nextflow/nextflow.config \
                run nextflow/clinvarbitration.nf \
                --ref_fa {ref_fa} \
                --output_dir {job.output} \
                --gff3 {gff3}
            """,
        )

        # copy the outputs back, in one smooooooth motion
        get_batch().write_output(job.output, str(outputs['release.tar.gz']).removesuffix('.release.tar.gz'))

        return self.make_outputs(multicohort, data=outputs, jobs=job)
