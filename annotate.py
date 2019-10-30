"""
Example annotation pipeline.

1. Load contigs.

2. Run seqstats.

3. Call genes.

4. COGProfiling.

"""
import argparse
from os import path, listdir
from abacat import Genome, timer_wrapper, CONFIG
from dfa_lib_python import (
    Dataflow,
    Attribute,
    AttributeType,
    add_transformation,
    start_task,
    end_task,
)


@timer_wrapper
def main(dataflow_tag, input_dir, db, transformations=False, blast="blastn"):
    """
    :param dataflow_tag:
    :param contigs_file:
    :param db:
    :return:
    """

    dataflow_tag = dataflow_tag
    df = Dataflow(dataflow_tag)

    def create_schemas():
        """
        This should only be executed if the schema hasn't been created in the database.
        """
        tf0_input_attr = [Attribute("contigs_file", AttributeType.FILE)]
        tf0_output_attr = [
            Attribute("genome_instance", AttributeType.TEXT),
            Attribute("genome_directory", AttributeType.TEXT),
            Attribute("genome_name", AttributeType.TEXT),
        ]
        tf0, tf0_input, tf0_output = add_transformation(
            df, "GenomeObject", tf0_input_attr, tf0_output_attr
        )

        tf1_input_attr = [Attribute("contigs_file", AttributeType.FILE)]
        tf1_output_attr = [
            Attribute(label, AttributeType.NUMERIC)
            for label in (
                "total_n",
                "total_seq",
                "avg_seq",
                "median_seq",
                "n50",
                "Min_seq",
                "max_seq",
            )
        ]
        tf1, tf1_input, tf1_output = add_transformation(
            df, "SequenceStatistics", tf1_input_attr, tf1_output_attr, tf0_output
        )

        tf2_input_attr = [Attribute("contigs_file", AttributeType.FILE)]
        tf2_output_attr = [
            Attribute(label, AttributeType.FILE)
            for label in (
                "GeneCalling_genes",
                "GeneCalling_proteins",
                "GeneCalling_cds",
                "GeneCalling_scores",
            )
        ]
        tf2, tf2_input, tf2_output = add_transformation(
            df, "GeneCalling", tf2_input_attr, tf2_output_attr, tf0_output
        )

        tf3_input_attr = [
            Attribute("GeneCalling_genes", AttributeType.FILE),
            Attribute("Database", AttributeType.TEXT),
            Attribute("Database_Path", AttributeType.FILE),
        ]
        tf3_output_attr = [
            Attribute("Database", AttributeType.TEXT),
            Attribute("Blast_xml_out", AttributeType.FILE),
            Attribute("Blast_seqs", AttributeType.FILE),
            Attribute("No_of_genes", AttributeType.NUMERIC),
        ]
        tf3, tf3_input, tf3_output = add_transformation(
            df, "COGProfiling", tf3_input_attr, tf3_output_attr, tf2_output
        )


        # Save all TFs in the dataflow
        df.save()

    if transformations:
        create_schemas()

    ix = 0
    for contigs_file in listdir(input_dir):
        contigs_file = path.join(input_dir, contigs_file)
        if contigs_file.endswith("_genomic.fna"):
            # Task 0 - GenomeObject
            t0, t0_input = start_task(ix, dataflow_tag, "GenomeObject", [contigs_file])
            genome = Genome(contigs_file)
            t0_output = end_task(
                t0, "GenomeObject", [genome.__repr__(), genome.directory, genome.name]
            )
            ix += 1

            # Task 1 - SequenceStatistics
            t1, t1_input = start_task(
                ix, dataflow_tag, "SequenceStatistics", [genome.files["contigs"]]
            )
            genome.load_seqstats()
            genome.seqstats
            t1_output = end_task(
                t1, "SequenceStatistics", [int(i) for i in genome.seqstats.values()]
            )
            ix += 1

            # Task 2 - GeneCalling
            t2, t2_input = start_task(
                ix, dataflow_tag, "GeneCalling", [genome.files["contigs"]]
            )
            genome.run_prodigal()
            t2_output = end_task(
                t2, "GeneCalling", list(genome.files["prodigal"].values())
            )
            ix += 1

            # Task 3 - COGProfiling
            t3, t3_input = start_task(
                ix,
                dataflow_tag,
                "COGProfiling",
                [genome.files["prodigal"]["genes"], db, CONFIG["db"][db]],
            )
            genome.blast_seqs(db=db, blast=blast)
            t3_output = end_task(
                t3,
                "COGProfiling",
                [
                    db,
                    genome.files[db]["xml"],
                    genome.files[db]["annotation"],
                    len(genome.geneset[db]["records"]),
                ],
            )
            ix += 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Annotation pipeline. Starts with a contig file."
    )
    parser.add_argument(
        "-df",
        "--dataflow_tag",
        type=str,
        help="Name of your dataflow. Must be short and with no spaces or special characters.",
        default="bactools",
    )
    parser.add_argument(
        "-i",
        "--input",
        help="Input file or directory. Must be valid FASTA contigs files (post-assembly).",
    )
    parser.add_argument(
        "-db",
        "--database",
        type=str,
        help="Database name. Must be in bactools.CONFIG.py db parameter.",
    )
    parser.add_argument(
        "-b", "--blast", type=str, help="Blast type to be used.", default="blastn"
    )
    parser.add_argument(
        "-t",
        "--transformation",
        type=bool,
        default=False,
        help="Whether or not to create the database schemas. This should only be true if it is the first time running the script.",
    )
    args = parser.parse_args()
    main(args.dataflow_tag, args.input, args.database, args.transformation, args.blast)
