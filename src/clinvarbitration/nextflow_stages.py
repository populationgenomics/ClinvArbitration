from argparse import ArgumentParser
from typing import TYPE_CHECKING

from cpg_flow.stage import MultiCohortStage, stage
from cpg_flow.workflow import run_workflow

from clinvarbitration.jobs.ClinvArbitrationNextflow import clinvarbitration_nextflow
from clinvarbitration.stages import get_output_folder, populate_job_meta

if TYPE_CHECKING:
    from cpg_flow.stage import StageInput, StageOutput
    from cpg_flow.targets.multicohort import MultiCohort
    from cpg_utils import Path


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
        outputs = self.expected_outputs(multicohort)
        job = clinvarbitration_nextflow(output_root=str(outputs['release.tar.gz']).removesuffix('.release.tar.gz'))
        return self.make_outputs(multicohort, data=outputs, jobs=job)


def main(dry_run: bool = False):
    run_workflow(
        stages=[
            ClinvarbitrationNextflow,
        ],
        dry_run=dry_run,
    )


def cli_main():
    parser = ArgumentParser()
    parser.add_argument('--dry-run', action='store_true', help='Print the commands that would be run')
    args = parser.parse_args()
    main(dry_run=args.dry_run)


if __name__ == '__main__':
    cli_main()
