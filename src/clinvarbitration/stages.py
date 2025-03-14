from datetime import datetime
from functools import cache
from os.path import join

from cpg_utils import Path, to_path
from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import get_batch

from cpg_flow.stage import MultiCohortStage, StageInput, StageOutput, stage
from cpg_flow.targets import MultiCohort


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


@stage
class CopyLatestClinvarFiles(MultiCohortStage):
    def expected_outputs(self, mc: MultiCohort) -> dict[str, Path]:
        return {
            'submission_file': get_output_folder() / 'submission_summary.txt.gz',
            'variant_file': get_output_folder() / 'variant_summary.txt.gz',
        }

    def queue_jobs(self, mc: MultiCohort, inputs: StageInput) -> StageOutput:
        """
        run a wget copy of the relevant files into GCP
        """
        bash_job = get_batch().new_bash_job('CopyLatestClinvarFiles')
        bash_job.image(config_retrieve(['workflow', 'driver_image']))

        directory = 'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/'
        sub_file = 'submission_summary.txt.gz'
        var_file = 'variant_summary.txt.gz'

        bash_job.command(f'wget -q {directory}{sub_file} -O {bash_job.subs}')
        bash_job.command(f'wget -q {directory}{var_file} -O {bash_job.vars}')

        outputs = self.expected_outputs(mc)

        get_batch().write_output(bash_job.subs, str(outputs['submission_file']))
        get_batch().write_output(bash_job.vars, str(outputs['variant_file']))

        return self.make_outputs(data=outputs, jobs=bash_job, target=mc)


@stage(required_stages=[CopyLatestClinvarFiles], analysis_type='clinvarbitration', analysis_keys=['clinvar_decisions'])
class GenerateNewClinvarSummary(MultiCohortStage):
    def expected_outputs(self, mc: MultiCohort) -> dict[str, Path]:
        """
        a couple of files and a HT as Paths
        """
        return {
            'clinvar_decisions': get_output_folder() / 'clinvar_decisions.ht.tar.gz',
            'snv_vcf': get_output_folder() / 'clinvar_decisions.vcf.bgz',
        }

    def queue_jobs(self, mc: MultiCohort, inputs: StageInput) -> StageOutput:
        # relatively RAM intensive, short running task
        job = get_batch().new_job('GenerateNewClinvarSummary')
        job.image(config_retrieve(['workflow', 'driver_image'])).memory('highmem').cpu('2')

        # get the expected outputs
        outputs = self.expected_outputs(mc)

        if sites_to_blacklist := config_retrieve(['workflow', 'site_blacklist'], []):
            blacklist_sites = ' '.join(f'"{site}"' for site in sites_to_blacklist)
            blacklist_string = f' -b {blacklist_sites}'
        else:
            blacklist_string = ''

        var_file = get_batch().read_input(inputs.as_str(mc, CopyLatestClinvarFiles, 'variant_file'))
        sub_file = get_batch().read_input(inputs.as_str(mc, CopyLatestClinvarFiles, 'submission_file'))

        job.declare_resource_group(
            output={
                'ht.tar.gz': '{root}.ht.tar.gz',
                'vcf.bgz': '{root}.vcf.bgz',
                'vcf.bgz.tbi': '{root}.vcf.bgz.tbi',
            },
        )

        # resummary is an entrypoint alias for scripts/resummarise_clinvar.py
        job.command(f'resummary -v {var_file} -s {sub_file} -o {job.output} --minimal {blacklist_string}')

        # don't tar from current location, we'll catch all the tmp pathing
        job.command(f'mv {job.output}.ht clinvar_decisions.ht && tar -czf {job.output}.ht.tar.gz clinvar_decisions.ht')

        # selectively copy back some outputs
        get_batch().write_output(job.output, str(outputs['clinvar_decisions']).removesuffix('.ht.tar.gz'))

        return self.make_outputs(target=mc, data=outputs, jobs=job)


@stage(required_stages=[GenerateNewClinvarSummary])
class AnnotateClinvarSnvsWithBcftools(MultiCohortStage):
    """
    take the vcf output from the clinvar stage, and apply consequence annotations
    this uses BCFtools for speed and re-deployability, instead of VEP
    """

    def expected_outputs(self, mc: MultiCohort) -> Path:
        return get_output_folder() / 'clinvar_decisions.annotated.tsv'

    def queue_jobs(self, mc: MultiCohort, inputs: StageInput) -> StageOutput:
        outputs = self.expected_outputs(mc)

        # need a genome reference - mandatory argument
        ref_fa = get_batch().read_input(config_retrieve(['workflow', 'ref_fa']))
        snv_vcf = get_batch().read_input(inputs.as_str(mc, GenerateNewClinvarSummary, 'snv_vcf'))

        genome_build = config_retrieve(['workflow', 'genome_build'], 'GRCh38')
        if genome_build not in ['GRCh37', 'GRCh38']:
            raise ValueError('Unsupported genome build specified')

        gff3 = f'/clinvarbitration/bcftools_data/{genome_build}.gff3.gz'

        job = get_batch().new_job('AnnotateClinvarSnvsWithBcftools', attributes={'tool': 'bcftools'})
        job.image(config_retrieve(['workflow', 'driver_image']))
        job.storage('10G')

        # -g is the GFF3 file, -f is the reference fasta
        # --local-csq is required to apply non-phase aware annotation
        # --force is required to use annotations without phase data
        job.command(
            f"""bcftools csq --force --local-csq -f {ref_fa} -g {gff3} {snv_vcf} -o out.vcf""",
        )

        # split the bcftools CSQ fields, filter to missense, and write out a tab-delimited file
        # -d - duplicate, writes each transcript onto a new line
        # -s :missense - only keep Consequence==missense variants
        # -f - format string - tab delimited, Transcript, Amino Acid Change, ClinVar allele ID, ClinVar gold stars
        job.command(
            f'bcftools +split-vep out.vcf -d -s :missense -f "%transcript\t%amino_acid_change\t%allele_id\t%gold_stars\n" > {job.output}',  # noqa: E501
        )

        get_batch().write_output(job.output, str(outputs))

        return self.make_outputs(target=mc, jobs=job, data=outputs)


@stage(required_stages=[AnnotateClinvarSnvsWithBcftools])
class Pm5TableGeneration(MultiCohortStage):
    def expected_outputs(self, mc: MultiCohort) -> dict[str, Path]:
        """
        a single HT, compressed into a tarball
        """
        return {
            'pm5_ht': get_output_folder() / 'clinvar_decisions.pm5.ht.tar.gz',
            'pm5_json': get_output_folder() / 'clinvar_decisions.pm5.json',
        }

    def queue_jobs(self, mc: MultiCohort, inputs: StageInput) -> StageOutput:
        job = get_batch().new_job('Pm5TableGeneration')
        job.image(config_retrieve(['workflow', 'driver_image']))

        # get the expected outputs
        outputs = self.expected_outputs(mc)

        annotated_snvs = get_batch().read_input(inputs.as_str(mc, AnnotateClinvarSnvsWithBcftools))

        job.declare_resource_group(output={'pm5.ht.tar.gz': '{root}.pm5.ht.tar.gz', 'json': '{root}.json'})

        # write both HT and JSON outputs to the same root location
        job.command(f'pm5_table -i {annotated_snvs} -o {job.output}')

        # compress the HT and remove as a single file
        job.command(
            f'mv {job.output}.ht clinvar_decisions.pm5.ht && '
            f'tar -czf {job.output}.pm5.ht.tar.gz clinvar_decisions.pm5.ht',
        )

        # write both outputs together
        get_batch().write_output(job.output, str(outputs['pm5_json']).removesuffix('.pm5.json'))

        return self.make_outputs(target=mc, data=outputs, jobs=job)


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

    def expected_outputs(self, multicohort: MultiCohort) -> Path:
        return get_output_folder() / 'clinvar_decisions.release.tar.gz'

    def queue_jobs(self, multicohort: MultiCohort, inputs: StageInput) -> StageOutput:
        """
        Localise all the previously generated data into a folder
        tarball it, and write out as a single file
        """
        tar_output = self.expected_outputs(multicohort)

        job = get_batch().new_job('PackageForRelease')
        job.storage('10G')
        job.image(config_retrieve(['workflow', 'driver_image']))

        # find paths to the previous outputs
        pm5 = inputs.as_dict(multicohort, Pm5TableGeneration)
        clinvar_decisions = inputs.as_str(multicohort, GenerateNewClinvarSummary, 'clinvar_decisions')

        annotated_tsv = get_batch().read_input(inputs.as_str(multicohort, AnnotateClinvarSnvsWithBcftools))
        pm5_ht = get_batch().read_input(str(pm5['pm5_ht']))
        pm5_json = get_batch().read_input(str(pm5['pm5_json']))
        decisions_ht = get_batch().read_input(clinvar_decisions)

        # create a new folder, move the files into it
        # unpack all the already compressed HTs into it
        # the compress once, containing all files for distribution
        job.command(
            f"""
            mkdir clinvarbitration_data
            mv {annotated_tsv} clinvarbitration_data/clinvar_decisions.annotated.tsv
            mv {pm5_json} clinvarbitration_data/clinvar_decisions.pm5.json
            tar -xzf {pm5_ht} -C clinvarbitration_data
            tar -xzf {decisions_ht} -C clinvarbitration_data
            tar -czf {job.output} \
                clinvarbitration_data/clinvar_decisions.annotated.tsv \
                clinvarbitration_data/clinvar_decisions.pm5.json \
                clinvarbitration_data/clinvar_decisions.ht \
                clinvarbitration_data/clinvar_decisions.pm5.ht
        """,
        )
        get_batch().write_output(job.output, str(tar_output))
        return self.make_outputs(multicohort, data=tar_output, jobs=job)
