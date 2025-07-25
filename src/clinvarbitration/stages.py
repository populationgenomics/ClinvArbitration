from datetime import datetime
from functools import cache
from os.path import join
from typing import TYPE_CHECKING

from cpg_utils import Path, to_path
from cpg_utils.config import config_retrieve

from cpg_flow.stage import MultiCohortStage, stage

from clinvarbitration import __version__ as clinvarbitration_version
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
    sits in a bucket accessible to the operating dataset, in a folder named "clinvarbitration/YY-MM"
    """

    return to_path(
        join(
            config_retrieve(['storage', 'default', 'default']),
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

    print(f'Generating output meta for {output_file}')
    return {'image': config_retrieve(['workflow', 'driver_image']), 'version': clinvarbitration_version}


@stage
class CopyLatestClinvarFiles(MultiCohortStage):
    """
    This stage runs a wget on two separate resource files:
    1. ClinVar submissions (every individual clinvar submission, submitter, rating...)
    2. ClinVar variants (every individual clinvar variant, position, gene, etc...)

    These are localised in preparation for the next stage, which will generate a re-summary of the clinvar data
    """

    def expected_outputs(self, mc: 'MultiCohort') -> 'dict[str, Path]':
        return {
            'submission_file': get_output_folder() / 'submission_summary.txt.gz',
            'variant_file': get_output_folder() / 'variant_summary.txt.gz',
        }

    def queue_jobs(self, mc: 'MultiCohort', inputs: 'StageInput') -> 'StageOutput':
        outputs = self.expected_outputs(mc)

        bash_job = copy_latest_files(
            submissions=outputs['submission_file'],
            variants=outputs['variant_file'],
        )

        return self.make_outputs(data=outputs, jobs=bash_job, target=mc)


@stage(required_stages=[CopyLatestClinvarFiles])
class GenerateNewClinvarSummary(MultiCohortStage):
    """
    Use the downloaded ClinVar raw data, and re-assess all submissions at each locus
    This is done by grouping all submissions by their ClinVar Allele ID,
    filtering out old/irrelevant submissions,
    then re-assessing the remaining submissions to create a more decisive classification
    """

    def expected_outputs(self, mc: 'MultiCohort') -> 'dict[str, Path]':
        return {
            'clinvar_decisions': get_output_folder() / 'clinvar_decisions.ht.tar.gz',
            'snv_vcf': get_output_folder() / 'clinvar_decisions.vcf.bgz',
        }

    def queue_jobs(self, mc: 'MultiCohort', inputs: 'StageInput') -> 'StageOutput':
        outputs = self.expected_outputs(mc)

        # identify input files from CopyLatestClinvarFiles (raw ClinVar data)
        var_file = inputs.as_str(mc, CopyLatestClinvarFiles, 'variant_file')
        sub_file = inputs.as_str(mc, CopyLatestClinvarFiles, 'submission_file')

        job = generate_new_summary(
            var_file=var_file,
            sub_file=sub_file,
            output_root=str(outputs['clinvar_decisions']).removesuffix('.ht.tar.gz'),
        )
        return self.make_outputs(target=mc, data=outputs, jobs=job)


@stage(required_stages=[GenerateNewClinvarSummary])
class AnnotateClinvarSnvsWithBcftools(MultiCohortStage):
    """
    Take the vcf output from the clinvar stage, and apply consequence annotations
    this uses BCFtools for speed and re-deployability, instead of VEP
    """

    def expected_outputs(self, mc: 'MultiCohort') -> 'Path':
        return get_output_folder() / 'clinvar_decisions.annotated.tsv'

    def queue_jobs(self, mc: 'MultiCohort', inputs: 'StageInput') -> 'StageOutput':
        output = self.expected_outputs(mc)

        snv_vcf = inputs.as_str(mc, GenerateNewClinvarSummary, 'snv_vcf')

        job = annotate_clinvar_snvs(
            snv_vcf=snv_vcf,
            output=output,
        )
        return self.make_outputs(target=mc, jobs=job, data=output)


@stage(required_stages=[AnnotateClinvarSnvsWithBcftools])
class Pm5TableGeneration(MultiCohortStage):
    """
    Reads in the annotated variant data (in TSV format), and generates a PM5 table
    For each missense variant, we collect all other pathogenic missense variants affecting the same codon
    This is output as an HT, and a JSON file of the raw representation
    """

    def expected_outputs(self, mc: 'MultiCohort') -> 'dict[str, Path]':
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
    Takes the data created so far, and packages it up for release
    This includes the re-summarised data, the annotated variants, and the PM5 table/JSON representation
    They are exported as a single tarball, which should be uploaded to the release page monthly
    """

    def expected_outputs(self, multicohort: 'MultiCohort') -> 'Path':
        return get_output_folder() / 'clinvar_decisions.release.tar.gz'

    def queue_jobs(self, multicohort: 'MultiCohort', inputs: 'StageInput') -> 'StageOutput':
        """
        Localise all the previously generated data into a folder
        tarball it, and write out as a single file
        """
        output = self.expected_outputs(multicohort)

        annotated_variants_tsv = inputs.as_path(multicohort, AnnotateClinvarSnvsWithBcftools)
        clinvar_decisions = inputs.as_str(multicohort, GenerateNewClinvarSummary, 'clinvar_decisions')
        pm5 = inputs.as_dict(multicohort, Pm5TableGeneration)

        job = package_data_for_release(
            annotated_tsv=annotated_variants_tsv,
            pm5_json=pm5['pm5_json'],
            pm5_ht=pm5['pm5_ht'],
            clinvar_decisions=clinvar_decisions,
            output=output,
        )

        return self.make_outputs(multicohort, data=output, jobs=job)
