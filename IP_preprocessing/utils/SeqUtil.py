import copy
import logging
import os
import re
import tempfile
from multiprocessing import Process
from pathlib import Path
import numpy as np
import socket

from Bio import AlignIO
from Bio import Phylo
from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Phylo.BaseTree import Clade, Tree
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from typing import List, Dict, Tuple


from .Enum import ChainType
from .FileUtil import detach_format, TableData, TableHeader, read_file_fast, write_file
from .ProcessUtil import run_igblast, run_igblast_new
from .ProcessUtil import run_IMPre

# New import for blast_run
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbiblastpCommandline
#from Bio.Blast import NCBIStandalone

from itertools import product

#_protein_alphabet = IUPACProtein()
#_protein_letters = str(_protein_alphabet.letters)


hostname = socket.gethostname()


def is_nucleotide_or_protein(sequence: str):
    """
    return TRUE if the input sequence is a dna sequence,
    FALSE if it is a protein sequence,
    None otherwise.
    :param sequence:
    :return:
    """
    if all(ch.upper() in "AGTC-. " for ch in sequence):
        return True
    elif all(ch.upper() in _protein_letters + "-. X" for ch in sequence):
        return False
    else:
        return None


class PTMFinder(object):
    """Post Translation Modification sites finder."""

    def __init__(self):
        """"""
        self.re_deamidation = re.compile("(NS|NG|NH)")
        self.re_aspartate_isomerization = re.compile("(DG|DS|DT|DD|DA)")
        self.re_N_glycosylation = re.compile("(N[^P][ST])")
        self.re_cleavage = re.compile("(DP|DQ|NS)")
        self.re_oxidation = re.compile("([MW])")
        self.re_free_cysteine = re.compile("(C)")
        self.re_sulfation = re.compile("(Y)")
        self.re_methylation = re.compile("([KR])")

    def find_diamidation(self, sequence) -> List:
        """
        Find deamidation sites from the given sequence.

        * Deamidation Sites
        Asparagine residues as potential deamidation sites were analysed in the context of
        both the amino and carboxy adjacent amino acids using a method based on
        Robinson & Robinson, PNAS (2001) 98,3, p944-949.
        No potential deamidation sites with T1/2 <25 days were identified within the VL.

        :param sequence: Target amino acid sequence to find.
        :return: List of match objects that refer the sites found.
        """

        if type(sequence) is str:
            return [_ for _ in self.re_deamidation.finditer(sequence)]
        elif type(sequence) is Seq and type(sequence.alphabet) in (
            ExtendedIUPACProtein,
            IUPACProtein,
        ):
            return [_ for _ in self.re_deamidation.finditer(str(sequence))]
        else:
            return []

    def find_aspartate_isomerization(self, sequence) -> List:
        """
        Find aspartate isomerization sites from the given sequence.

        * Aspartate Isomerization Sites
        Aspartate isomerization sites were predicted by analysing the sequence
        for known isomerization motifs: DG, DS, DT, DD or DA.
        Two potential isomerization sites were identified in the VL sequence.
        Aspartate 91 and Aspartate 95A, both of which are located in CDR3.

        :param sequence: Target amino acid sequence to find.
        :return: List of match objects that refer the sites found.
        """

        if type(sequence) is str:
            return [_ for _ in self.re_aspartate_isomerization.finditer(sequence)]
        elif type(sequence) is Seq and type(sequence.alphabet) in (
            ExtendedIUPACProtein,
            IUPACProtein,
        ):
            return [_ for _ in self.re_aspartate_isomerization.finditer(str(sequence))]
        else:
            return []

    def find_n_glycosylation(self, sequence) -> List:
        """
        Find N glycosylation sites from the given sequence.

        * N Glycosylation Sites
        Sequences were analysed based on the consensus N-linked glycosylation motif:
        - N-X-S/T - (where X can be any amino acid except Proline)
        No potential N-linked glycosylation sites were identified within the VL domain.

        :param sequence: Target amino acid sequence to find.
        :return: List of match objects that refer the sites found.
        """

        if type(sequence) is str:
            return [_ for _ in self.re_N_glycosylation.finditer(sequence)]
        elif type(sequence) is Seq and type(sequence.alphabet) in (
            ExtendedIUPACProtein,
            IUPACProtein,
        ):
            return [_ for _ in self.re_N_glycosylation.finditer(str(sequence))]
        else:
            return []

    def find_cleavage(self, sequence) -> List:
        """
        Find cleavage sites from the given sequence.

        :param sequence: Target amino acid sequence to find.
        :return: List of match objects that refer the sites found.
        """

        if type(sequence) is str:
            return [_ for _ in self.re_cleavage.finditer(sequence)]
        elif type(sequence) is Seq and type(sequence.alphabet) in (
            ExtendedIUPACProtein,
            IUPACProtein,
        ):
            return [_ for _ in self.re_cleavage.finditer(str(sequence))]
        else:
            return []

    def find_oxidation(self, sequence) -> List:
        """
        Find oxidation sites from the given sequence.

        * Oxidation Sites
        Structural models were prepared and methionine and tryptophan residues were identified and assessed
        to determine whether they are likely to be surface exposed (and therefore candidates for oxidation).

        :param sequence: Target amino acid sequence to find.
        :return: List of match objects that refer the sites found.
        """

        if type(sequence) is str:
            return [_ for _ in self.re_oxidation.finditer(sequence)]
        elif type(sequence) is Seq and type(sequence.alphabet) in (
            ExtendedIUPACProtein,
            IUPACProtein,
        ):
            return [_ for _ in self.re_oxidation.finditer(str(sequence))]
        else:
            return []

    def find_free_cysteine(self, sequence) -> List:
        """
        Find free cysteine sites from the given sequence.

        * Free Cysteine Sites Sequences were analysed to identify unpaired cysteine residues. VL domains normally
        contain two cysteines which form a disulphide bond in the folded molecule. Additional cysteines would be
        expected to be detrimental to folding and potentially cause issues such as aggregation. No unpaired cysteines
        were identified within the VL domain.

        :param sequence: Target amino acid sequence to find.
        :return: List of match objects that refer the sites found.
        """

        if type(sequence) is str:
            return [_ for _ in self.re_free_cysteine.finditer(sequence)]
        elif type(sequence) is Seq and type(sequence.alphabet) in (
            ExtendedIUPACProtein,
            IUPACProtein,
        ):
            return [_ for _ in self.re_free_cysteine.finditer(str(sequence))]
        else:
            return []

    def find_sulfation(self, sequence) -> List:
        """
        Find sulfation sites from the given sequence.

        * Sulfation Sites Tyrosine sulfation sites were predicted using the method of Monigatti F. et al.,
        Bioinformatics 18:769-770 (2002). No sulfated tyrosine sites were predicted

        :param sequence: Target amino acid sequence to find.
        :return: List of match objects that refer the sites found.
        """

        if type(sequence) is str:
            return [_ for _ in self.re_sulfation.finditer(sequence)]
        elif type(sequence) is Seq and type(sequence.alphabet) in (
            ExtendedIUPACProtein,
            IUPACProtein,
        ):
            return [_ for _ in self.re_sulfation.finditer(str(sequence))]
        else:
            return []

    def find_methylation(self, sequence) -> List:
        """
        Find methylation sites from the given sequence.

        * Methylation Sites
        Prediction of methylation of lysine and arginine residues was performed using the method of
        Chen et al., Nucleic Acids Research  34 (Issue suppl 2): 249-253 (2006)
        One methylated lysine site was predicted – Lysine 66 (located in Framework 3)
        PSRFSGS - K - SGSTHTL
        No methylated arginine sites were predicted.

        :param sequence: Target amino acid sequence to find.
        :return: List of match objects that refer the sites found.
        """

        if type(sequence) is str:
            return [_ for _ in self.re_methylation.finditer(sequence)]
        elif type(sequence) is Seq and type(sequence.alphabet) in (
            ExtendedIUPACProtein,
            IUPACProtein,
        ):
            return [_ for _ in self.re_methylation.finditer(str(sequence))]
        else:
            return []

    @staticmethod
    def find_custom(sequence, re_compiled) -> List:
        """
        Find custom subsequence  sites from the given sequence.

        :param sequence: Target amino acid sequence to find.
        :param re_compiled:
        :return: List of match objects that refer the sites found.
        """

        if type(sequence) is str:
            return [_ for _ in re_compiled.finditer(sequence)]
        elif type(sequence) is Seq and type(sequence.alphabet) in (
            ExtendedIUPACProtein,
            IUPACProtein,
        ):
            return [_ for _ in re_compiled.finditer(str(sequence))]
        else:
            return []


class Helper(object):
    def __init__(self):
        """"""
        self.reverse_sequence = {
            "A": "T",
            "G": "C",
            "T": "A",
            "C": "G",
            "W": "W",
            "S": "S",
            "M": "K",
            "K": "M",
            "R": "Y",
            "Y": "R",
            "B": "V",
            "V": "B",
            "D": "H",
            "H": "D",
            "N": "N",
            " ": " ",
            ".": ".",
        }
        self.standard_table = CodonTable.unambiguous_dna_by_id[1]
        self.standard_codon_map = copy.deepcopy(self.standard_table.forward_table)
        for triplet in self.standard_table.stop_codons:
            self.standard_codon_map[triplet] = "*"
        self.standard_codon_map["   "] = " "
        self.standard_codon_map["..."] = "."
        self.dna_symbol_map = {
            "A": "A",
            "C": "C",
            "G": "G",
            "T": "T",
            "U": "U",
            "W": "[AT]",
            "S": "[CG]",
            "M": "[AC]",
            "K": "[GT]",
            "R": "[AG]",
            "Y": "[CT]",
            "B": "[CGT]",
            "D": "[AGT]",
            "H": "[ACT]",
            "V": "[ACG]",
            "N": "[ACGT]",
        }

        self.IUPAC_code = {
            "A": "A",
            "G": "G",
            "C": "C",
            "T": "T",
            "U": "U",
            "W": ["A", "T"],
            "S": ["C", "G"],
            "M": ["A", "C"],
            "K": ["G", "T"],
            "R": ["A", "G"],
            "Y": ["C", "T"],
            "B": ["C", "G", "T"],
            "D": ["A", "G", "T"],
            "H": ["A", "C", "T"],
            "V": ["A", "C", "G"],
            "N": ["A", "C", "G", "T"],
        }

        # amber codon TAG->Q
        self.amber_codon_map = copy.deepcopy(self.standard_table.forward_table)
        for triplet in self.standard_table.stop_codons:
            self.amber_codon_map[triplet] = "*"
        self.amber_codon_map["TAG"] = "Q"

    def extend_degenerate_dna(self, sequence: str):
        """return list of all possible sequences given an degenerate DNA input"""
        # r = []
        # for i in product(*[self.IUPAC_code[j] for j in sequence]):
        #     r.append("".join(i))
        r = [
            "".join(i) for i in product(*[self.IUPAC_code[j] for j in sequence.upper()])
        ]
        return r

    def to_reverse_complement(self, sequence: str) -> str:
        """"""
        return "".join(self.reverse_sequence[s] for s in sequence)[::-1]

    def to_regex_string_from_dna(self, sequence: str, is_forward: bool = True) -> str:
        """
        Symbol Table
        A	Adenine	A
        C	Cytosine		C
        G	Guanine			G
        T	Thymine				T
        U	Uracil				U
        W	Weak	A			T
        S	Strong		C	G
        M	aMino	A	C
        K	Keto			G	T
        R	puRine	A		G
        Y	pYrimidine		C		T
        B	not A (B comes after A)		C	G	T
        D	not C (D comes after C)	A		G	T
        H	not G (H comes after G)	A	C		T
        V	not T (V comes after T and U)	A	C	G
        N	any Nucleotide (not a gap)	A	C	G	T
        :param sequence:
        :return:
        """

        if is_forward:
            return "".join([self.dna_symbol_map[c] for c in sequence.upper()])
        else:
            return "".join(
                [
                    self.dna_symbol_map[c]
                    for c in self.to_reverse_complement(sequence.upper())
                ]
            )

    @staticmethod
    def hamming_distance(s1: str, s2: str) -> int:
        """Calculate the Hamming distance between two strings"""
        if len(s1) != len(s2):
            return max(len(s1), len(s2))

        return sum((c1 != c2 and c1 != "." and c2 != ".") for c1, c2 in zip(s1, s2))

    @staticmethod
    def get_diff(ref: str, sub: str):
        if len(ref) != len(sub):
            return None

        return [(i, c2) for i, (c1, c2) in enumerate(zip(ref, sub)) if c1 != c2]

    @staticmethod
    def hamming_distance_from_diff(df1: str, df2: str):
        if len(df1) == 0:
            return len(df2)
        elif len(df2) == 0:
            return len(df1)

        i1, i2, n_same_position, n_same_pair = 0, 0, 0, 0
        while i1 < len(df1) and i2 < len(df2):
            if df1[i1][0] > df2[i2][0]:
                i2 += 1
            elif df1[i1][0] < df2[i2][0]:
                i1 += 1
            else:
                n_same_position += 1
                if df1[i1][1] == df2[i2][1]:
                    n_same_pair += 1
                i1 += 1
                i2 += 1

        return len(df1) + len(df2) - n_same_position - n_same_pair

    @staticmethod
    def levenshtein_distance(s1: str, s2: str) -> int:
        return editdistance.eval(s1, s2)

    @staticmethod
    def to_bitstream(dna_sequence: str):
        stream1 = 0x00
        stream2 = 0x00
        for nucl in dna_sequence:
            if nucl == "A":
                pass
            elif nucl == "G":
                stream1 += 1
            elif nucl == "C":
                stream2 += 1
            elif nucl == "T":
                stream1 += 1
                stream2 += 1
            else:
                return None
            stream1 = stream1 << 1
            stream2 = stream2 << 1
        return stream1, stream2

    @staticmethod
    def hamming_distance_from_bitstream(stream_pair1, stream_pair2):
        return bin(
            (stream_pair1[0] ^ stream_pair2[0]) | (stream_pair1[1] ^ stream_pair2[1])
        ).count("1")

    def translate(self, dna: str, codon_map: Dict[str, str] = None, **kwargs) -> str:
        if not codon_map:
            codon_map = self.standard_codon_map

        if codon_map == "amber":
            codon_map = self.amber_codon_map

        aa = ""
        for i in range(0, len(dna), 3):
            triplet = dna[i : i + 3]
            try:
                aa += codon_map[triplet]
            except KeyError:
                if len(triplet) < 3:
                    aa += "X"
                else:
                    aa += "X"
            except Exception as e:
                print(e)

        return aa

    @staticmethod
    def RepresentsInt(s):
        try:
            int(s)
            return True
        except ValueError:
            return False

    def annotate_regions_old(
        self,
        chain_type: ChainType,
        sequence_list: List[str],
        seq_type: str = "Ig",
        **kwargs
    ) -> TableData:
        """
        :param chain_type:
        :param sequence_list:
        :return:
        """

        # check sequence type
        is_nucl = is_nucleotide_or_protein(sequence_list[0])
        if is_nucl == None or is_nucl == False:
            raise RuntimeError(
                "parameter sequence_list is not a nucleotide sequence list"
            )

        if chain_type in [
            ChainType.HUMAN_HEAVY,
            ChainType.HUMAN_LIGHT,
            ChainType.CHICKEN_HEAVY,
            ChainType.CHICKEN_LIGHT,
            ChainType.RAT,
            ChainType.RABBIT,
            ChainType.MOUSE_BALBC,
            ChainType.MOUSE_C57BL6,
        ]:
            tp_query = tempfile.NamedTemporaryFile("wt", suffix=".fasta", delete=False)
            if "file_blast" in kwargs:
                file_blast = kwargs["file_blast"]
            else:
                tp = tempfile.NamedTemporaryFile("wt")
                file_blast = tp.name
                tp.close()

            if "domain_system" in kwargs:
                domain_system = kwargs["domain_system"]
            else:
                domain_system = "kabat"

            if "blast_kwargs" in kwargs:
                blast_kwargs = kwargs["blast_kwargs"]
            else:
                blast_kwargs = {"num_threads": 10}

            if "seq_id_list" in kwargs:
                seq_id_list = kwargs["seq_id_list"]
                for i, (seq, seq_id) in enumerate(zip(sequence_list, seq_id_list)):
                    tp_query.write(">%s\n" % seq_id)
                    tp_query.write(seq + "\n")
            else:
                for i, seq in enumerate(sequence_list):
                    tp_query.write(">seq_%d\n" % i)
                    tp_query.write(seq + "\n")

            tp_query.close()
            run_igblast(
                chain_type,
                tp_query.name,
                file_blast,
                seq_type,
                domain_system,
                **blast_kwargs
            )
            print("[LOG] make %s done" % file_blast)

            result_list = self.parse_igblast(
                chain_type=chain_type,
                igblast_output=file_blast,
                domain_system=domain_system,
            )

            for sequence, blast_result in zip(sequence_list, result_list):
                aa_seq = self.translate(sequence)
                if len(sequence) % 3 == 0:
                    frame = "In"
                else:
                    frame = "Out"
                blast_result["query"] = sequence
                blast_result["frame"] = frame
                blast_result["stop codon"] = "*" in aa_seq
                if "FR4 from" in blast_result:
                    blast_result["FR4"] = sequence[
                        blast_result["FR4 from"] : blast_result["FR4 from"] + 33
                    ]
                    blast_result["FR4 to"] = blast_result["FR4 from"] + 33
                else:
                    blast_result["FR4"] = "N/A"
                if chain_type in [
                    ChainType.HUMAN_HEAVY,
                    ChainType.HUMAN_LIGHT,
                    ChainType.RABBIT,
                ]:
                    for r in ["FR1", "FR2", "FR3", "CDR1", "CDR2"]:
                        if ("%s from" % r) in blast_result and (
                            "%s to" % r
                        ) in blast_result:
                            blast_result[r] = sequence[
                                blast_result["%s from" % r] : blast_result["%s to" % r]
                            ]
                        else:
                            blast_result[r] = "N/A"
                elif chain_type in [ChainType.RAT]:
                    for r in ["FR1", "FR2", "FR3", "CDR1", "CDR2", "CDR3"]:
                        if ("%s from" % r) in blast_result and (
                            "%s to" % r
                        ) in blast_result:
                            blast_result[r] = sequence[
                                blast_result["%s from" % r] : blast_result["%s to" % r]
                            ]
                        else:
                            blast_result[r] = "N/A"

            return result_list
        else:
            raise RuntimeError("chain type error")

    def annotate_regions(
        self,
        chain_type: ChainType,
        sequence_list: List[str],
        annotation_list: List[str],
        **kwargs
    ) -> TableData:
        """
        :param chain_type:
        :param sequence_list:
        :return:
        """

        # check sequence type
        is_nucl = is_nucleotide_or_protein(sequence_list[0])

        if is_nucl == None or is_nucl == False:
            raise RuntimeError(
                "parameter sequence_list is not a nucleotide sequence list"
            )

        if chain_type in [
            ChainType.HUMAN_HEAVY,
            ChainType.HUMAN_LIGHT,
            ChainType.CHICKEN_HEAVY,
            ChainType.CHICKEN_LIGHT,
            ChainType.RAT_HEAVY,
            ChainType.RAT_LIGHT,
            ChainType.RABBIT_HEAVY,
            ChainType.RABBIT_KAPPA,
            ChainType.MOUSE_BALBC_HEAVY,
            ChainType.MOUSE_BALBC_LIGHT,
            ChainType.MOUSE_C57BL6_HEAVY,
            ChainType.MOUSE_C57BL6_LIGHT,
            ChainType.MOUSE_HEAVY,
            ChainType.MOUSE_LIGHT,
            ChainType.ALPACA_HEAVY,
        ]:
            seq_type = "Ig"
            # domain_system = 'kabat'
            domain_system = "imgt"

        elif chain_type in [
            ChainType.HUMAN_BETA,
            ChainType.HUMAN_ALPHA,
            ChainType.HUMAN_DELTA,
            ChainType.HUMAN_GAMMA,
            ChainType.MOUSE_BETA,
            ChainType.MOUSE_ALPHA,
            ChainType.MOUSE_GAMMA,
            ChainType.MOUSE_DELTA,
        ]:
            seq_type = "TCR"
            domain_system = "imgt"

        else:
            raise RuntimeError("chain type error")

        if "file_blast" in kwargs:
            file_blast = kwargs["file_blast"]
        else:
            tp = tempfile.NamedTemporaryFile("wt")
            file_blast = tp.name
            tp.close()

        if "blast_kwargs" in kwargs:
            blast_kwargs = kwargs["blast_kwargs"]
        else:
            blast_kwargs = {}

        tp_query = tempfile.NamedTemporaryFile("wt", suffix=".fasta", delete=False)
        if "seq_id_list" in kwargs:
            seq_id_list = kwargs["seq_id_list"]
            for i, (seq, seq_id) in enumerate(zip(sequence_list, seq_id_list)):
                tp_query.write(">%s\n" % seq_id)
                tp_query.write(seq + "\n")
        else:
            for i, seq in enumerate(sequence_list):
                tp_query.write(">seq_%d\n" % i)
                tp_query.write(seq + "\n")
        tp_query.close()

        # run_igblast(chain_type, tp_query.name, file_blast, seq_type, domain_system, **blast_kwargs)
        run_igblast_new(
            chain_type,
            tp_query.name,
            file_blast,
            seq_type,
            domain_system,
            **blast_kwargs
        )
        logging.info("[LOG] make %s done" % file_blast)

        result_list = self.parse_igblast(
            chain_type=chain_type,
            igblast_output=file_blast,
            domain_system=domain_system,
        )

        ### combine annotation information with input_data_list (isotyped input file)
        for sequence, igblast_parsed_data in zip(sequence_list, result_list):
            tmp_dict = {"sequence": sequence}
            igblast_parsed_data.update(tmp_dict)

            ## check not annotated cases and fill it with default info (N/A or blank)
            if "sequence_aa" not in igblast_parsed_data:
                igblast_parsed_data["sequence_aa"] = ""
            if igblast_parsed_data["productive"] == "Yes":
                igblast_parsed_data["productive"] = "T"
            else:
                igblast_parsed_data["productive"] = "F"

            # add cdr1, 2
            if igblast_parsed_data["rev_comp"] == "T":
                sequence = self.to_reverse_complement(sequence)
            if "fr4_from" in igblast_parsed_data:
                igblast_parsed_data["fr4"] = sequence[
                    igblast_parsed_data["fr4_from"] : igblast_parsed_data["fr4_from"]
                    + 33
                ]
            for r in ["fr1", "fr2", "fr3", "cdr1", "cdr2"]:
                if ("%s_from" % r) in igblast_parsed_data and (
                    "%s_to" % r
                ) in igblast_parsed_data:
                    igblast_parsed_data[r] = sequence[
                        igblast_parsed_data["%s_from" % r] : igblast_parsed_data[
                            "%s_to" % r
                        ]
                    ]
                    igblast_parsed_data[r + "_aa"] = self.translate(
                        igblast_parsed_data[r]
                    )
                else:
                    igblast_parsed_data[r] = "N/A"
                    igblast_parsed_data[r + "_aa"] = "X"

            # add annotation_list data
            for annotation in annotation_list:
                if annotation not in igblast_parsed_data:
                    igblast_parsed_data[annotation] = "N/A"
                else:
                    continue

        return result_list

    def annotate_regions_chicken(
        self, seq_list, v_type, result_dir="/home/IP-team/data/blast_db/blast_temp/"
    ):
        query_dir = os.path.join(result_dir, "query")
        out_dir = os.path.join(result_dir, "out")

        if os.path.exists(query_dir) is not True:
            os.makedirs(query_dir)

        if os.path.exists(out_dir) is not True:
            os.makedirs(out_dir)

        v_type = v_type.lower()

        with open(os.path.join(query_dir, "query.fasta"), "w") as fasta_writer:
            for i, e_seq in enumerate(seq_list):
                fasta_writer.write(">" + str(i) + "\n")
                fasta_writer.write(e_seq + "\n")

        db_dir = "/home/IP-team/data/blast_db/fr_db_chicken_white_leghorn"

        region_list = ["fr1", "fr2", "fr3", "fr4"]

        # Variables used in region extraction & Initialization
        ext_index = []
        cdr_list = []

        for i in range(len(seq_list)):
            temp = {
                "fr1_e": None,
                "fr2_s": None,
                "fr2_e": None,
                "fr3_s": None,
                "fr3_e": None,
                "fr4_s": None,
            }
            ext_index.append(temp)

        # Parameters for blast_run
        E_VAL = 5
        NUM_ALIGN = 10
        THREADS = 20
        FORMAT_EXE = "/Tools/ncbi-blast-2.7.1+/bin/blastn"

        for region in region_list:
            db_name = "_".join([v_type, region])
            db_file = os.path.join(db_dir, db_name)
            out_name = "query_" + db_name + "_blasted.txt"
            out_file = os.path.join(out_dir, out_name)

            print("Blast run for " + db_name + " is started.")

            # Note that word_size = 5
            cline = NcbiblastnCommandline(
                num_threads=THREADS,
                query=os.path.join(query_dir, "query.fasta"),
                db=db_file,
                evalue=0.1**E_VAL,
                out=out_file,
                gapopen=1,
                gapextend=2,
                word_size=5,
                num_descriptions=NUM_ALIGN,
                num_alignments=NUM_ALIGN,
            )

            os.system(FORMAT_EXE + str(cline)[len("blastn") :])

            print("Blast run for " + db_name + " is completed.\n")

            handle = open(out_file, "r")
            blast_parser = NCBIStandalone.BlastParser()
            iterator = NCBIStandalone.Iterator(handle, blast_parser)

            for i, each in enumerate(iterator):
                if len(each.alignments) == 0:
                    continue

                for alignment in each.alignments:
                    hsps = alignment.hsps[0]

                    s_point = hsps.query_start
                    e_point = hsps.query_end

                    if region == "fr1":
                        if (
                            ext_index[i]["fr1_e"] == None
                            and hsps.sbjct_end == alignment.length
                        ):
                            ext_index[i]["fr1_e"] = e_point
                            break
                    elif region == "fr2":
                        if ext_index[i]["fr2_s"] == None and hsps.sbjct_start == 1:
                            ext_index[i]["fr2_s"] = s_point - 1
                        if (
                            ext_index[i]["fr2_e"] == None
                            and hsps.sbjct_end == alignment.length
                        ):
                            ext_index[i]["fr2_e"] = e_point
                        if (
                            ext_index[i]["fr2_s"] != None
                            and ext_index[i]["fr2_e"] != None
                        ):
                            break
                    elif region == "fr3":
                        if ext_index[i]["fr3_s"] == None and hsps.sbjct_start == 1:
                            ext_index[i]["fr3_s"] = s_point - 1
                        if (
                            ext_index[i]["fr3_e"] == None
                            and hsps.sbjct_end == alignment.length
                        ):
                            ext_index[i]["fr3_e"] = e_point
                        if (
                            ext_index[i]["fr3_s"] != None
                            and ext_index[i]["fr3_e"] != None
                        ):
                            break
                    elif region == "fr4":
                        if ext_index[i]["fr4_s"] == None and hsps.sbjct_start == 1:
                            ext_index[i]["fr4_s"] = s_point - 1
                            break

        for seq, ext in zip(seq_list, ext_index):
            if ext["fr1_e"] == None or ext["fr2_s"] == None:
                cdr1 = "E/F"
            else:
                cdr1 = seq[ext["fr1_e"] : ext["fr2_s"]]

            if ext["fr2_e"] == None or ext["fr3_s"] == None:
                cdr2 = "E/F"
            else:
                cdr2 = seq[ext["fr2_e"] : ext["fr3_s"]]

            if ext["fr3_e"] == None or ext["fr4_s"] == None:
                cdr3 = "E/F"
            else:
                cdr3 = seq[ext["fr3_e"] : ext["fr4_s"]]

            cdr_list.append({"cdr1": cdr1, "cdr2": cdr2, "cdr3": cdr3})

        fail_list = []
        for i, e_cdr in enumerate(cdr_list):
            if e_cdr["cdr3"] == "E/F":
                fail_list.append(i)

        return cdr_list, fail_list

    def annotate_regions_chicken_loose(
        self,
        seq_list,
        v_type,
        e_val=5,
        result_dir="/home/IP-team/data/blast_db/blast_temp/",
    ):
        query_dir = os.path.join(result_dir, "query")
        out_dir = os.path.join(result_dir, "out")

        if os.path.exists(query_dir) is not True:
            os.makedirs(query_dir)

        if os.path.exists(out_dir) is not True:
            os.makedirs(out_dir)

        v_type = v_type.lower()

        with open(os.path.join(query_dir, "query.fasta"), "w") as fasta_writer:
            for i, e_seq in enumerate(seq_list):
                fasta_writer.write(">" + str(i) + "\n")
                fasta_writer.write(e_seq + "\n")

        db_dir = "/home/IP-team/data/blast_db/fr_db_chicken_white_leghorn"

        region_list = ["fr1", "fr2", "fr3", "fr4"]

        # Variables used in region extraction & Initialization
        ext_index = []
        cdr_list = []

        for i in range(len(seq_list)):
            temp = {
                "fr1_e": None,
                "fr2_s": None,
                "fr2_e": None,
                "fr3_s": None,
                "fr3_e": None,
                "fr4_s": None,
            }
            ext_index.append(temp)

        # Parameters for blast_run
        E_VAL = e_val
        NUM_ALIGN = 1
        THREADS = 20
        FORMAT_EXE = "/Tools/ncbi-blast-2.7.1+/bin/blastn"

        for region in region_list:
            db_name = "_".join([v_type, region])
            db_file = os.path.join(db_dir, db_name)
            out_name = "query_" + db_name + "_blasted.txt"
            out_file = os.path.join(out_dir, out_name)

            print("Blast run for " + db_name + " is started.")

            # Note that word_size = 10
            cline = NcbiblastnCommandline(
                num_threads=THREADS,
                query=os.path.join(query_dir, "query.fasta"),
                db=db_file,
                evalue=0.1**E_VAL,
                out=out_file,
                gapopen=1,
                gapextend=2,
                word_size=10,
                num_descriptions=NUM_ALIGN,
                num_alignments=NUM_ALIGN,
            )

            os.system(FORMAT_EXE + str(cline)[len("blastn") :])

            print("Blast run for " + db_name + " is completed.\n")

            handle = open(out_file, "r")
            blast_parser = NCBIStandalone.BlastParser()
            iterator = NCBIStandalone.Iterator(handle, blast_parser)

            for i, each in enumerate(iterator):
                if len(each.alignments) == 0:
                    continue

                for alignment in each.alignments:
                    hsps = alignment.hsps[0]

                    s_point = hsps.query_start - hsps.sbjct_start + 1
                    e_point = hsps.query_end + alignment.length - hsps.sbjct_end

                    if region == "fr1":
                        ext_index[i]["fr1_e"] = e_point
                        break
                    elif region == "fr2":
                        ext_index[i]["fr2_s"] = s_point - 1
                        ext_index[i]["fr2_e"] = e_point
                        break
                    elif region == "fr3":
                        ext_index[i]["fr3_s"] = s_point - 1
                        ext_index[i]["fr3_e"] = e_point
                        break
                    elif region == "fr4":
                        ext_index[i]["fr4_s"] = s_point - 1
                        break

        for seq, ext in zip(seq_list, ext_index):
            if ext["fr1_e"] == None or ext["fr2_s"] == None:
                cdr1 = "E/F"
            else:
                cdr1 = seq[ext["fr1_e"] : ext["fr2_s"]]

            if ext["fr2_e"] == None or ext["fr3_s"] == None:
                cdr2 = "E/F"
            else:
                cdr2 = seq[ext["fr2_e"] : ext["fr3_s"]]

            if ext["fr3_e"] == None or ext["fr4_s"] == None:
                cdr3 = "E/F"
            else:
                cdr3 = seq[ext["fr3_e"] : ext["fr4_s"]]

            cdr_list.append({"cdr1": cdr1, "cdr2": cdr2, "cdr3": cdr3})

        return cdr_list

    def parse_igblast(self, chain_type: ChainType, igblast_output: str, **kwargs):
        sHelper = Helper()
        result_list = []
        # result = None

        if "domain_system" in kwargs:
            domain_system = kwargs["domain_system"]
        else:
            # domain_system = 'kabat'
            domain_system = "imgt"

        ### parse igblast result txt file.
        igblast_parsed_dict = {}
        with open(igblast_output, "r") as igb_handle:
            # for the case that domain system is "imgt"
            if domain_system == "imgt":
                # include BCR & TCR
                if chain_type in [
                    ChainType.HUMAN_HEAVY,
                    ChainType.HUMAN_LIGHT,
                    ChainType.RAT_HEAVY,
                    ChainType.RAT_LIGHT,
                    ChainType.RABBIT_HEAVY,
                    ChainType.RABBIT_KAPPA,
                    ChainType.MOUSE_BALBC_HEAVY,
                    ChainType.MOUSE_BALBC_LIGHT,
                    ChainType.MOUSE_C57BL6_HEAVY,
                    ChainType.MOUSE_C57BL6_LIGHT,
                    ChainType.ALPACA_HEAVY,
                    ChainType.MOUSE_HEAVY,
                    ChainType.MOUSE_LIGHT,
                    ChainType.HUMAN_BETA,
                    ChainType.HUMAN_ALPHA,
                    ChainType.HUMAN_DELTA,
                    ChainType.HUMAN_GAMMA,
                    ChainType.MOUSE_BETA,
                    ChainType.MOUSE_ALPHA,
                    ChainType.MOUSE_DELTA,
                    ChainType.MOUSE_GAMMA,
                ]:
                    # for the case that domain system is "imgt"
                    while True:
                        line = igb_handle.readline()
                        if not line:
                            break

                        if line[: len("Query= ")] == "Query= ":
                            query_id = line[len("Query= ") :].strip()

                            junction = ""
                            junction_aa = ""
                            junction_start = 0
                            junction_end = 0
                            CDR3_NT = ""
                            sequence_alignment = ""
                            seq_aa = ""

                            result = {
                                "query": line[len("Query= ") : -1],
                                "hit": False,
                                "v_call": "N/A",
                                "j_call": "N/A",
                                "cdr3": "N/A",
                                "cdr3_aa": "X",
                                "junction": "N/A",
                                "junction_aa": "X",
                                "rev_comp": "F",
                                "productive": "F",
                                "v_cigar": "",
                                "j_cigar": "",
                                "v_alignment_length": "N/A",
                                "v_alignment_mutation": "N/A",
                                "sequence_alignment": "N/A",
                                "germline_alignment": "N/A",
                            }
                            if chain_type in [
                                ChainType.HUMAN_HEAVY,
                                ChainType.RABBIT_HEAVY,
                                ChainType.MOUSE_C57BL6_HEAVY,
                                ChainType.MOUSE_BALBC_HEAVY,
                                ChainType.MOUSE_HEAVY,
                                ChainType.HUMAN_BETA,
                                ChainType.HUMAN_DELTA,
                                ChainType.MOUSE_BETA,
                                ChainType.MOUSE_DELTA,
                            ]:
                                result["d_call"] = "N/A"
                                result["d_cigar"] = ""

                            while True:
                                line = igb_handle.readline()

                                if "No hits found" in line:
                                    result["hit"] = False
                                    igblast_parsed_dict[query_id] = result
                                    break

                                elif (
                                    "Note that your query represents the minus strand"
                                    in line
                                ):
                                    result["rev_comp"] = "T"

                                elif (
                                    line[: len("Sub-region sequence details")]
                                    == "Sub-region sequence details"
                                ):
                                    next_line_sp = igb_handle.readline().split("\t")
                                    if next_line_sp[0] == "CDR3":
                                        result["hit"] = True
                                        CDR3_NT = next_line_sp[1]
                                        result["cdr3"] = CDR3_NT
                                        result["cdr3_aa"] = sHelper.translate(
                                            result["cdr3"]
                                        )
                                        result["cdr3_start"] = int(next_line_sp[3])
                                        junction_start = int(next_line_sp[3]) - 4
                                        if junction_start != 0:
                                            result["junction_start"] = junction_start
                                        result["cdr3_end"] = int(next_line_sp[4])
                                        junction_end = int(next_line_sp[4]) + 3
                                        if junction_end != 0:
                                            result["junction_end"] = junction_end
                                        result["fr4_from"] = int(next_line_sp[4])

                                elif line[: len("Alignments")] == "Alignments":
                                    next_line_sp = igb_handle.readline()
                                    next_line_sp = igb_handle.readline()
                                    next_line_sp = igb_handle.readline().split(" ")
                                    next_line_trimmed = []
                                    for l in range(len(next_line_sp)):
                                        if next_line_sp[l] != "":
                                            next_line_trimmed.append(next_line_sp[l])
                                    query_line = []
                                    query_start_where = []
                                    query_finish_where = []
                                    V_germline_mismatch = []
                                    D_germline_mismatch = []
                                    J_germline_mismatch = []
                                    V_mis_ratio = []
                                    D_mis_ratio = []
                                    J_mis_ratio = []
                                    v_germ_start_where = []
                                    v_germ_finish_where = []
                                    d_germ_start_where = []
                                    d_germ_finish_where = []
                                    j_germ_start_where = []
                                    j_germ_finish_where = []
                                    sequence_a = []
                                    dog = True
                                    while dog:
                                        next_line_trimmed = []
                                        before_line_sp = next_line_sp
                                        next_line_sp = igb_handle.readline().split(" ")
                                        for l in range(len(next_line_sp)):
                                            if next_line_sp[l] != "":
                                                next_line_trimmed.append(
                                                    next_line_sp[l]
                                                )

                                        before_line_trimmed = []
                                        for l in range(len(before_line_sp)):
                                            if before_line_sp[l] != "":
                                                before_line_trimmed.append(
                                                    before_line_sp[l]
                                                )

                                        findQu = re.compile("Query_")
                                        findQue = findQu.search(next_line_trimmed[0])

                                        findper = re.compile("%")
                                        if len(next_line_trimmed) >= 2:
                                            findperc = findper.search(
                                                next_line_trimmed[1]
                                            )

                                        if findQue:
                                            query_start_where.append(
                                                next_line_trimmed[1]
                                            )
                                            query_line.append(next_line_trimmed[2])
                                            query_finish_where.append(
                                                next_line_trimmed[3].replace("\n", "")
                                            )
                                            seq_a = ""
                                            for l in range(len(before_line_trimmed)):
                                                seq_a = seq_a + before_line_trimmed[l]
                                            sequence_a.append(seq_a.replace("\n", ""))

                                        elif next_line_trimmed[0] == "V" and findperc:
                                            if len(next_line_trimmed) >= 7:
                                                if len(V_mis_ratio) == 0:
                                                    V_mis_ratio.append(
                                                        next_line_trimmed[1:3]
                                                    )
                                                v_germ_start_where.append(
                                                    next_line_trimmed[4]
                                                )
                                                V_germline_mismatch.append(
                                                    next_line_trimmed[5]
                                                )
                                                v_germ_finish_where.append(
                                                    next_line_trimmed[6].replace(
                                                        "\n", ""
                                                    )
                                                )

                                        elif next_line_trimmed[0] == "D" and findperc:
                                            if len(next_line_trimmed) >= 7:
                                                if len(D_mis_ratio) == 0:
                                                    D_mis_ratio.append(
                                                        next_line_trimmed[1:3]
                                                    )
                                                d_germ_start_where.append(
                                                    next_line_trimmed[4]
                                                )
                                                D_germline_mismatch.append(
                                                    next_line_trimmed[5]
                                                )
                                                d_germ_finish_where.append(
                                                    next_line_trimmed[6].replace(
                                                        "\n", ""
                                                    )
                                                )

                                        elif next_line_trimmed[0] == "J" and findperc:
                                            if len(next_line_trimmed) >= 7:
                                                if len(J_mis_ratio) == 0:
                                                    J_mis_ratio.append(
                                                        next_line_trimmed[1:3]
                                                    )
                                                j_germ_start_where.append(
                                                    next_line_trimmed[4]
                                                )
                                                J_germline_mismatch.append(
                                                    next_line_trimmed[5]
                                                )
                                                j_germ_finish_where.append(
                                                    next_line_trimmed[6].replace(
                                                        "\n", ""
                                                    )
                                                )

                                        elif next_line_trimmed[0] == "Lambda":
                                            dog = False
                                    seq_aa = ""
                                    for k in range(len(sequence_a)):
                                        seq_aa = seq_aa + sequence_a[k]
                                    result["sequence_aa"] = seq_aa
                                    result["query_line"] = query_line
                                    for a in range(len(query_line)):
                                        sequence_alignment = (
                                            sequence_alignment + query_line[a]
                                        )
                                    result["sequence_alignment"] = sequence_alignment
                                    seq_without = ""
                                    if len(sequence_alignment) != 0:
                                        for z in range(len(sequence_alignment)):
                                            if sequence_alignment[z] != "-":
                                                seq_without = (
                                                    seq_without + sequence_alignment[z]
                                                )
                                    if len(seq_without) != 0 and junction_start != 0:
                                        junction = seq_without[
                                            junction_start
                                            - int(query_start_where[0])
                                            + 1 : junction_end
                                            - int(query_start_where[0])
                                            + 1
                                        ]
                                    if CDR3_NT != "":
                                        result["junction"] = junction
                                        result["junction_aa"] = sHelper.translate(
                                            junction
                                        )
                                    else:
                                        result["junction"] = ""
                                        result["junction_aa"] = ""
                                    v_finish_from_zero = 0
                                    if len(v_germ_start_where) != 0:
                                        v_finish_from_zero = int(
                                            v_germ_finish_where[
                                                len(v_germ_finish_where) - 1
                                            ]
                                        ) - int(v_germ_start_where[0])
                                    v_cigar_mismatch = []
                                    v_cigar_mismatches = []
                                    j_cigar_mismatches = []
                                    d_cigar_mismatches = []
                                    if (
                                        len(sequence_alignment)
                                        >= v_finish_from_zero + 1
                                    ) and v_finish_from_zero != 0:
                                        for a in range(v_finish_from_zero + 1):
                                            if sequence_alignment[a] == "-":
                                                v_cigar_mismatch.append([a, "0"])
                                    if len(v_germ_start_where) != 0:
                                        for a in range(len(v_germ_start_where)):
                                            for i in range(
                                                int(v_germ_finish_where[a])
                                                - int(v_germ_start_where[a])
                                                + 1
                                            ):
                                                if V_germline_mismatch[a][i] != ".":
                                                    v_cigar_mismatch.append(
                                                        [
                                                            i
                                                            + int(v_germ_start_where[a])
                                                            - int(
                                                                v_germ_start_where[0]
                                                            ),
                                                            V_germline_mismatch[a][i],
                                                        ]
                                                    )
                                    v_cigar_mismatch.sort(key=lambda x: x[0])
                                    for kk in range(len(v_cigar_mismatch)):
                                        isit = 0
                                        for ll in range(len(v_cigar_mismatches)):
                                            if (
                                                v_cigar_mismatch[kk][0]
                                                == v_cigar_mismatches[ll][0]
                                            ):
                                                isit += 1
                                        if isit == 0:
                                            v_cigar_mismatches.append(
                                                v_cigar_mismatch[kk]
                                            )

                                    d_finish_from_zero = 0
                                    if len(d_germ_start_where) != 0:
                                        d_finish_from_zero = int(
                                            d_germ_finish_where[
                                                len(d_germ_finish_where) - 1
                                            ]
                                        ) - int(d_germ_start_where[0])

                                    d_cigar_mismatch = []
                                    d_line_pos = []
                                    if (len(d_germ_start_where)) != 0:
                                        for k in range(90):
                                            if D_germline_mismatch[0][k] != "-":
                                                d_line_po = k
                                                d_line_pos.append(k)
                                                break

                                    seq_before_d_num = 0
                                    if (
                                        d_finish_from_zero != 0
                                        and len(v_germ_finish_where) != 0
                                        and len(sequence_alignment) != 0
                                    ):
                                        for q in range(len(v_germ_start_where) - 1):
                                            # seq_before_d_num : 90, 180, ...
                                            seq_before_d_num = seq_before_d_num + 90
                                        for a in range(len(d_cigar_mismatch)):
                                            if (
                                                sequence_alignment[
                                                    seq_before_d_num
                                                    + d_line_pos[0]
                                                    + 90 * a
                                                ]
                                                == "-"
                                            ):
                                                d_cigar_mismatch.append(
                                                    [
                                                        seq_before_d_num
                                                        + d_line_pos[0]
                                                        + 90 * a,
                                                        "0",
                                                    ]
                                                )
                                    if len(d_germ_start_where) == 1:
                                        for i in range(
                                            int(d_germ_finish_where[0])
                                            - int(d_germ_start_where[0])
                                        ):
                                            if (
                                                D_germline_mismatch[0][
                                                    d_line_pos[0] + i
                                                ]
                                                != "."
                                            ):
                                                d_cigar_mismatch.append(
                                                    [
                                                        i
                                                        + seq_before_d_num
                                                        + d_line_pos[0],
                                                        D_germline_mismatch[0][
                                                            d_line_pos[0] + i
                                                        ],
                                                    ]
                                                )
                                    elif len(d_germ_start_where) >= 2:
                                        for a in range(90 - d_line_pos[0]):
                                            if (
                                                D_germline_mismatch[0][
                                                    d_line_pos[0] + a
                                                ]
                                                != "."
                                            ):
                                                d_cigar_mismatch.append(
                                                    [
                                                        seq_before_d_num
                                                        + d_line_pos[0]
                                                        + a,
                                                        D_germline_mismatch[0][
                                                            d_line_pos[0] + a
                                                        ],
                                                    ]
                                                )
                                        for a in range(
                                            int(d_germ_finish_where[1])
                                            - int(d_germ_start_where[0])
                                            + 1
                                            - 90
                                            + d_line_pos[0]
                                        ):
                                            if D_germline_mismatch[1][a] != ".":
                                                d_cigar_mismatch.append(
                                                    [
                                                        seq_before_d_num + 90 + a,
                                                        D_germline_mismatch[1][a],
                                                    ]
                                                )
                                    d_cigar_mismatch.sort(key=lambda x: x[0])
                                    for kk in range(len(d_cigar_mismatch)):
                                        isit = 0
                                        for ll in range(len(d_cigar_mismatches)):
                                            if (
                                                d_cigar_mismatch[kk][0]
                                                == d_cigar_mismatches[ll][0]
                                            ):
                                                isit += 1
                                        if isit == 0:
                                            d_cigar_mismatches.append(
                                                d_cigar_mismatch[kk]
                                            )

                                    j_finish_from_zero = 0
                                    if len(j_germ_start_where) != 0:
                                        j_finish_from_zero = int(
                                            j_germ_finish_where[
                                                len(j_germ_finish_where) - 1
                                            ]
                                        ) - int(j_germ_start_where[0])

                                    j_cigar_mismatch = []
                                    j_line_pos = []
                                    if (len(j_germ_start_where)) != 0:
                                        for k in range(90):
                                            if J_germline_mismatch[0][k] != "-":
                                                j_line_po = k
                                                j_line_pos.append(j_line_po)
                                                break

                                    seq_before_j_num = 0
                                    if (
                                        j_finish_from_zero != 0
                                        and len(v_germ_finish_where) != 0
                                    ):
                                        if (len(d_germ_finish_where)) <= 1:
                                            for q in range(len(v_germ_start_where) - 1):
                                                # seq_before_j_num : 90, 180, ...
                                                seq_before_j_num = seq_before_j_num + 90
                                        elif (len(d_germ_finish_where)) == 2:
                                            for q in range(len(v_germ_start_where)):
                                                # seq_before_j_num : 90, 180, ...
                                                seq_before_j_num = seq_before_j_num + 90
                                        for a in range(len(j_cigar_mismatch)):
                                            if (
                                                sequence_alignment[
                                                    seq_before_j_num
                                                    + j_line_pos[0]
                                                    + 90 * a
                                                ]
                                                == "-"
                                            ):
                                                j_cigar_mismatch.append(
                                                    [
                                                        seq_before_j_num
                                                        + j_line_pos[0]
                                                        + 90 * a,
                                                        "0",
                                                    ]
                                                )
                                    if len(j_germ_start_where) == 1:
                                        for i in range(
                                            int(j_germ_finish_where[0])
                                            - int(j_germ_start_where[0])
                                        ):
                                            if (
                                                J_germline_mismatch[0][
                                                    j_line_pos[0] + i
                                                ]
                                                != "."
                                            ):
                                                j_cigar_mismatch.append(
                                                    [
                                                        i
                                                        + seq_before_j_num
                                                        + j_line_pos[0],
                                                        J_germline_mismatch[0][
                                                            j_line_pos[0] + i
                                                        ],
                                                    ]
                                                )
                                    elif len(j_germ_start_where) >= 2:
                                        for a in range(90 - j_line_pos[0]):
                                            if (
                                                J_germline_mismatch[0][
                                                    j_line_pos[0] + a
                                                ]
                                                != "."
                                            ):
                                                j_cigar_mismatch.append(
                                                    [
                                                        seq_before_j_num
                                                        + j_line_pos[0]
                                                        + a,
                                                        J_germline_mismatch[0][
                                                            j_line_pos[0] + a
                                                        ],
                                                    ]
                                                )
                                        for pp in range(len(j_germ_start_where) - 1):
                                            for a in range(
                                                int(j_germ_finish_where[pp + 1])
                                                - int(j_germ_start_where[pp + 1])
                                                + 1
                                            ):
                                                if (
                                                    J_germline_mismatch[pp + 1][a]
                                                    != "."
                                                ):
                                                    j_cigar_mismatch.append(
                                                        [
                                                            seq_before_j_num
                                                            + 90 * (pp + 1)
                                                            + a,
                                                            J_germline_mismatch[1][a],
                                                        ]
                                                    )
                                    j_cigar_mismatch.sort(key=lambda x: x[0])

                                    for kk in range(len(j_cigar_mismatch)):
                                        isit = 0
                                        for ll in range(len(j_cigar_mismatches)):
                                            if (
                                                j_cigar_mismatch[kk][0]
                                                == j_cigar_mismatches[ll][0]
                                            ):
                                                isit += 1
                                        if isit == 0:
                                            j_cigar_mismatches.append(
                                                j_cigar_mismatch[kk]
                                            )

                                    germline_alignment = sequence_alignment
                                    germline_alignment_one = ""
                                    germline_alignment_two = ""
                                    germline_alignment_three = ""
                                    if (
                                        len(j_cigar_mismatches) != 0
                                        and len(germline_alignment) != 0
                                    ):
                                        for p in range(len(j_cigar_mismatches)):
                                            germline_alignment = (
                                                germline_alignment[
                                                    : j_cigar_mismatches[p][0]
                                                ]
                                                + j_cigar_mismatches[p][1]
                                                + germline_alignment[
                                                    j_cigar_mismatches[p][0] + 1 :
                                                ]
                                            )
                                    if (
                                        len(d_cigar_mismatches) != 0
                                        and len(germline_alignment) != 0
                                    ):
                                        for p in range(len(d_cigar_mismatches)):
                                            germline_alignment = (
                                                germline_alignment[
                                                    : d_cigar_mismatches[p][0]
                                                ]
                                                + d_cigar_mismatches[p][1]
                                                + germline_alignment[
                                                    d_cigar_mismatches[p][0] + 1 :
                                                ]
                                            )
                                    if (
                                        len(j_cigar_mismatch) != 0
                                        and len(germline_alignment) != 0
                                    ):
                                        for p in range(len(j_cigar_mismatch)):
                                            germline_alignment = (
                                                germline_alignment[
                                                    : j_cigar_mismatch[p][0]
                                                ]
                                                + j_cigar_mismatch[p][1]
                                                + germline_alignment[
                                                    j_cigar_mismatch[p][0] + 1 :
                                                ]
                                            )
                                    the_string = ""
                                    second_string = ""
                                    if (
                                        (len(d_germ_start_where)) == 0
                                        and (len(j_germ_finish_where)) != 0
                                        and len(v_germ_finish_where) != 0
                                    ):
                                        seq_before_j_num = 0
                                        if (
                                            j_finish_from_zero != 0
                                            and len(v_germ_finish_where) != 0
                                        ):
                                            if (len(d_germ_finish_where)) <= 1:
                                                for q in range(
                                                    len(v_germ_start_where) - 1
                                                ):
                                                    # seq_before_j_num : 90, 180, ...
                                                    seq_before_j_num = (
                                                        seq_before_j_num + 90
                                                    )
                                            elif (len(d_germ_finish_where)) == 2:
                                                for q in range(len(v_germ_start_where)):
                                                    # seq_before_j_num : 90, 180, ...
                                                    seq_before_j_num = (
                                                        seq_before_j_num + 90
                                                    )
                                        the_string = ""
                                        for u in range(
                                            seq_before_j_num
                                            + j_line_pos[0]
                                            - int(
                                                v_germ_finish_where[
                                                    (len(v_germ_finish_where) - 1)
                                                ]
                                            )
                                            + int(v_germ_start_where[0])
                                        ):
                                            the_string = the_string + "-"
                                        if len(germline_alignment) != 0:
                                            germline_alignment = (
                                                germline_alignment[
                                                    : int(
                                                        v_germ_finish_where[
                                                            len(v_germ_finish_where) - 1
                                                        ]
                                                    )
                                                    - int(v_germ_start_where[0])
                                                ]
                                                + the_string
                                                + germline_alignment[
                                                    seq_before_j_num + j_line_pos[0] :
                                                ]
                                            )
                                    elif (
                                        (len(d_germ_start_where)) != 0
                                        and len(j_germ_finish_where) != 0
                                        and len(v_germ_finish_where) != 0
                                    ):
                                        seq_before_d_num = 0
                                        for q in range(len(v_germ_start_where) - 1):
                                            # seq_before_d_num : 90, 180, ...
                                            seq_before_d_num = (
                                                seq_before_d_num
                                                + int(v_germ_finish_where[0])
                                                - int(v_germ_start_where[0])
                                                + 1
                                            )

                                        seq_before_j_num = 0
                                        if (
                                            j_finish_from_zero != 0
                                            and len(v_germ_finish_where) != 0
                                        ):
                                            if (len(d_germ_finish_where)) <= 1:
                                                for q in range(
                                                    len(v_germ_start_where) - 1
                                                ):
                                                    # seq_before_j_num : 90, 180, ...
                                                    seq_before_j_num = (
                                                        seq_before_j_num
                                                        + int(v_germ_finish_where[0])
                                                        - int(v_germ_start_where[0])
                                                        + 1
                                                    )
                                            elif (len(d_germ_finish_where)) == 2:
                                                for q in range(len(v_germ_start_where)):
                                                    # seq_before_j_num : 90, 180, ...
                                                    seq_before_j_num = (
                                                        seq_before_j_num
                                                        + int(v_germ_finish_where[0])
                                                        - int(v_germ_start_where[0])
                                                        + 1
                                                    )

                                        the_string = ""
                                        second_string = ""
                                        for u in range(
                                            seq_before_d_num
                                            + d_line_pos[0]
                                            - int(
                                                v_germ_finish_where[
                                                    (len(v_germ_finish_where) - 1)
                                                ]
                                            )
                                            + int(v_germ_start_where[0])
                                            - 1
                                        ):
                                            the_string = the_string + "-"
                                        for u in range(
                                            seq_before_j_num
                                            + j_line_pos[0]
                                            - seq_before_d_num
                                            - d_line_pos[0]
                                            - int(
                                                d_germ_finish_where[
                                                    len(d_germ_finish_where) - 1
                                                ]
                                            )
                                            + int(d_germ_start_where[0])
                                            - 1
                                        ):
                                            second_string = second_string + "-"
                                        germline_alignment_one = germline_alignment[
                                            : int(
                                                v_germ_finish_where[
                                                    len(v_germ_finish_where) - 1
                                                ]
                                            )
                                            - int(v_germ_start_where[0])
                                            + 1
                                        ]
                                        germline_alignment_two = germline_alignment[
                                            seq_before_d_num
                                            + d_line_pos[0] : seq_before_d_num
                                            + d_line_pos[0]
                                            + int(
                                                d_germ_finish_where[
                                                    len(d_germ_finish_where) - 1
                                                ]
                                            )
                                            - int(d_germ_start_where[0])
                                            + 1
                                        ]
                                        germline_alignment_three = germline_alignment[
                                            seq_before_j_num + j_line_pos[0] :
                                        ]
                                        germline_alignment = (
                                            germline_alignment_one
                                            + the_string
                                            + germline_alignment_two
                                            + second_string
                                            + germline_alignment_three
                                        )

                                    result["the_string"] = the_string
                                    result["second_string"] = second_string
                                    result["germline_alignment"] = germline_alignment
                                    result[
                                        "germline_alignment_1"
                                    ] = germline_alignment_one
                                    result[
                                        "germline_alignment_2"
                                    ] = germline_alignment_two
                                    result[
                                        "germline_alignment_3"
                                    ] = germline_alignment_three
                                    result["j_line_pos"] = j_line_pos

                                    v_cig = ""
                                    cur_pos_is = -1
                                    type_mis = 0
                                    v_cig_num = 0
                                    v_cigar = ""
                                    if len(v_cigar_mismatches) >= 1:
                                        cur_pos_is = -1
                                        for e in range(len(v_cigar_mismatches)):
                                            if v_cigar_mismatches[e][1] == "0":
                                                v_cig_num = (
                                                    int(v_cigar_mismatches[e][0])
                                                    - cur_pos_is
                                                    - 1
                                                )
                                                if v_cig_num == -1:
                                                    v_cig = v_cig
                                                elif v_cig_num != 0:
                                                    v_cig = (
                                                        v_cig
                                                        + str(v_cig_num)
                                                        + "="
                                                        + str(1)
                                                        + "D"
                                                    )
                                                    cur_pos_is = int(
                                                        v_cigar_mismatches[e][0]
                                                    )
                                                elif type_mis == 1:
                                                    v_cig = (
                                                        v_cig[: len(v_cig) - 2]
                                                        + str(
                                                            int(v_cig[len(v_cig) - 2])
                                                            + 1
                                                        )
                                                        + "D"
                                                    )
                                                    cur_pos_is = int(
                                                        v_cigar_mismatches[e][0]
                                                    )
                                                else:
                                                    v_cig = v_cig + str(1) + "D"
                                                    cur_pos_is = int(
                                                        v_cigar_mismatches[e][0]
                                                    )
                                                type_mis = 1
                                            elif v_cigar_mismatches[e][1] == "-":
                                                v_cig_num = (
                                                    int(v_cigar_mismatches[e][0])
                                                    - cur_pos_is
                                                    - 1
                                                )
                                                if v_cig_num == -1:
                                                    v_cig = v_cig
                                                elif v_cig_num != 0:
                                                    v_cig = (
                                                        v_cig
                                                        + str(v_cig_num)
                                                        + "="
                                                        + str(1)
                                                        + "I"
                                                    )
                                                    cur_pos_is = int(
                                                        v_cigar_mismatches[e][0]
                                                    )
                                                elif type_mis == 2:
                                                    v_cig = (
                                                        v_cig[: len(v_cig) - 2]
                                                        + str(
                                                            int(v_cig[len(v_cig) - 2])
                                                            + 1
                                                        )
                                                        + "I"
                                                    )
                                                    cur_pos_is = int(
                                                        v_cigar_mismatches[e][0]
                                                    )
                                                else:
                                                    v_cig = v_cig + str(1) + "I"
                                                    cur_pos_is = int(
                                                        v_cigar_mismatches[e][0]
                                                    )
                                                type_mis = 2
                                            else:
                                                v_cig_num = (
                                                    int(v_cigar_mismatches[e][0])
                                                    - cur_pos_is
                                                    - 1
                                                )
                                                if v_cig_num == -1:
                                                    v_cig = v_cig
                                                elif v_cig_num != 0:
                                                    v_cig = (
                                                        v_cig
                                                        + str(v_cig_num)
                                                        + "="
                                                        + str(1)
                                                        + "X"
                                                    )
                                                    cur_pos_is = int(
                                                        v_cigar_mismatches[e][0]
                                                    )
                                                elif type_mis == 3:
                                                    v_cig = (
                                                        v_cig[: len(v_cig) - 2]
                                                        + str(
                                                            int(v_cig[len(v_cig) - 2])
                                                            + 1
                                                        )
                                                        + "X"
                                                    )
                                                    cur_pos_is = int(
                                                        v_cigar_mismatches[e][0]
                                                    )
                                                else:
                                                    v_cig = v_cig + str(1) + "X"
                                                    cur_pos_is = int(
                                                        v_cigar_mismatches[e][0]
                                                    )
                                                type_mis = 3
                                    if len(v_germ_finish_where) != 0:
                                        if (
                                            cur_pos_is
                                            != int(
                                                v_germ_finish_where[
                                                    (len(v_germ_finish_where) - 1)
                                                ]
                                            )
                                            - int(v_germ_start_where[0])
                                            + 1
                                        ):
                                            v_cig_num = (
                                                int(
                                                    v_germ_finish_where[
                                                        (len(v_germ_finish_where) - 1)
                                                    ]
                                                )
                                                - int(v_germ_start_where[0])
                                                - cur_pos_is
                                            )
                                            v_cig = v_cig + str(v_cig_num) + "="
                                    v_cigar = v_cig

                                    d_cig = ""
                                    type_mis = 0
                                    d_cig_num = 0
                                    d_cigar = ""
                                    if (
                                        len(d_germ_finish_where) != 0
                                        and len(v_germ_finish_where) != 0
                                    ):
                                        cur_pos_is = (
                                            seq_before_d_num + d_line_pos[0] - 1
                                        )
                                        for e in range(len(d_cigar_mismatches)):
                                            if d_cigar_mismatches[e][1] == "0":
                                                d_cig_num = (
                                                    int(d_cigar_mismatches[e][0])
                                                    - cur_pos_is
                                                    - 1
                                                )
                                                if d_cig_num == -1:
                                                    d_cig = d_cig
                                                elif d_cig_num != 0:
                                                    d_cig = (
                                                        d_cig
                                                        + str(d_cig_num)
                                                        + "="
                                                        + str(1)
                                                        + "D"
                                                    )
                                                    cur_pos_is = int(
                                                        d_cigar_mismatches[e][0]
                                                    )
                                                elif type_mis == 1:
                                                    d_cig = (
                                                        d_cig[: len(d_cig) - 2]
                                                        + str(
                                                            int(d_cig[len(d_cig) - 2])
                                                            + 1
                                                        )
                                                        + "D"
                                                    )
                                                    cur_pos_is = int(
                                                        d_cigar_mismatches[e][0]
                                                    )
                                                else:
                                                    d_cig = d_cig + str(1) + "D"
                                                    cur_pos_is = int(
                                                        d_cigar_mismatches[e][0]
                                                    )
                                                type_mis = 1
                                            elif d_cigar_mismatches[e][1] == "-":
                                                d_cig_num = (
                                                    int(d_cigar_mismatches[e][0])
                                                    - cur_pos_is
                                                    - 1
                                                )
                                                if d_cig_num == -1:
                                                    d_cig = d_cig
                                                elif d_cig_num != 0:
                                                    d_cig = (
                                                        d_cig
                                                        + str(d_cig_num)
                                                        + "="
                                                        + str(1)
                                                        + "I"
                                                    )
                                                    cur_pos_is = int(
                                                        d_cigar_mismatches[e][0]
                                                    )
                                                elif type_mis == 2:
                                                    d_cig = (
                                                        d_cig[: len(d_cig) - 2]
                                                        + str(
                                                            int(d_cig[len(d_cig) - 2])
                                                            + 1
                                                        )
                                                        + "I"
                                                    )
                                                    cur_pos_is = int(
                                                        d_cigar_mismatches[e][0]
                                                    )
                                                else:
                                                    d_cig = d_cig + str(1) + "I"
                                                    cur_pos_is = int(
                                                        d_cigar_mismatches[e][0]
                                                    )
                                                type_mis = 2
                                            else:
                                                d_cig_num = (
                                                    int(d_cigar_mismatches[e][0])
                                                    - cur_pos_is
                                                    - 1
                                                )
                                                if d_cig_num == -1:
                                                    d_cig = d_cig
                                                elif d_cig_num != 0:
                                                    d_cig = (
                                                        d_cig
                                                        + str(d_cig_num)
                                                        + "="
                                                        + str(1)
                                                        + "X"
                                                    )
                                                    cur_pos_is = int(
                                                        d_cigar_mismatches[e][0]
                                                    )
                                                elif type_mis == 3:
                                                    d_cig = (
                                                        d_cig[: len(d_cig) - 2]
                                                        + str(
                                                            int(d_cig[len(d_cig) - 2])
                                                            + 1
                                                        )
                                                        + "X"
                                                    )
                                                    cur_pos_is = int(
                                                        d_cigar_mismatches[e][0]
                                                    )
                                                else:
                                                    d_cig = d_cig + str(1) + "X"
                                                    cur_pos_is = int(
                                                        d_cigar_mismatches[e][0]
                                                    )
                                                type_mis = 3
                                        if (
                                            cur_pos_is
                                            != int(
                                                d_germ_finish_where[
                                                    (len(d_germ_finish_where) - 1)
                                                ]
                                            )
                                            - int(d_germ_start_where[0])
                                            + seq_before_d_num
                                            + d_line_pos[0]
                                            + 1
                                        ):
                                            d_cig_num = (
                                                int(
                                                    d_germ_finish_where[
                                                        (len(d_germ_finish_where) - 1)
                                                    ]
                                                )
                                                - int(d_germ_start_where[0])
                                                - cur_pos_is
                                                + seq_before_d_num
                                                + d_line_pos[0]
                                            )
                                            d_cig = d_cig + str(d_cig_num) + "="
                                        d_cigar = d_cig
                                    elif len(d_germ_finish_where) == 0:
                                        d_cigar = ""

                                    j_cig = ""
                                    type_mis = 0
                                    j_cig_num = 0
                                    j_cigar = ""
                                    if (len(j_line_pos)) != 0:
                                        cur_pos_is = (
                                            seq_before_j_num + j_line_pos[0] - 1
                                        )
                                        for e in range(len(j_cigar_mismatches)):
                                            if j_cigar_mismatches[e][1] == "0":
                                                j_cig_num = (
                                                    int(j_cigar_mismatches[e][0])
                                                    - cur_pos_is
                                                    - 1
                                                )
                                                if j_cig_num == -1:
                                                    j_cig = j_cig
                                                elif j_cig_num != 0:
                                                    j_cig = (
                                                        j_cig
                                                        + str(j_cig_num)
                                                        + "="
                                                        + str(1)
                                                        + "D"
                                                    )
                                                    cur_pos_is = int(
                                                        j_cigar_mismatches[e][0]
                                                    )
                                                elif type_mis == 1:
                                                    j_cig = (
                                                        j_cig[: len(j_cig) - 2]
                                                        + str(
                                                            int(j_cig[len(j_cig) - 2])
                                                            + 1
                                                        )
                                                        + "D"
                                                    )
                                                    cur_pos_is = int(
                                                        j_cigar_mismatches[e][0]
                                                    )
                                                else:
                                                    j_cig = j_cig + str(1) + "D"
                                                    cur_pos_is = int(
                                                        j_cigar_mismatches[e][0]
                                                    )
                                                type_mis = 1
                                            elif j_cigar_mismatches[e][1] == "-":
                                                j_cig_num = (
                                                    int(j_cigar_mismatches[e][0])
                                                    - cur_pos_is
                                                    - 1
                                                )
                                                if j_cig_num == -1:
                                                    j_cig = j_cig
                                                elif j_cig_num != 0:
                                                    j_cig = (
                                                        j_cig
                                                        + str(j_cig_num)
                                                        + "="
                                                        + str(1)
                                                        + "I"
                                                    )
                                                    cur_pos_is = int(
                                                        j_cigar_mismatches[e][0]
                                                    )
                                                elif type_mis == 2:
                                                    j_cig = (
                                                        j_cig[: len(j_cig) - 2]
                                                        + str(
                                                            int(j_cig[len(j_cig) - 2])
                                                            + 1
                                                        )
                                                        + "I"
                                                    )
                                                    cur_pos_is = int(
                                                        j_cigar_mismatches[e][0]
                                                    )
                                                else:
                                                    j_cig = j_cig + str(1) + "I"
                                                    cur_pos_is = int(
                                                        j_cigar_mismatches[e][0]
                                                    )
                                                type_mis = 2
                                            else:
                                                j_cig_num = (
                                                    int(j_cigar_mismatches[e][0])
                                                    - cur_pos_is
                                                    - 1
                                                )
                                                if j_cig_num == -1:
                                                    j_cig = j_cig
                                                elif j_cig_num != 0:
                                                    j_cig = (
                                                        j_cig
                                                        + str(j_cig_num)
                                                        + "="
                                                        + str(1)
                                                        + "X"
                                                    )
                                                    cur_pos_is = int(
                                                        j_cigar_mismatches[e][0]
                                                    )
                                                elif type_mis == 3:
                                                    j_cig = (
                                                        j_cig[: len(j_cig) - 2]
                                                        + str(
                                                            int(j_cig[len(j_cig) - 2])
                                                            + 1
                                                        )
                                                        + "X"
                                                    )
                                                    cur_pos_is = int(
                                                        j_cigar_mismatches[e][0]
                                                    )
                                                else:
                                                    j_cig = j_cig + str(1) + "X"
                                                    cur_pos_is = int(
                                                        j_cigar_mismatches[e][0]
                                                    )
                                                type_mis = 3
                                        if (
                                            cur_pos_is
                                            != int(
                                                j_germ_finish_where[
                                                    (len(j_germ_finish_where) - 1)
                                                ]
                                            )
                                            - int(j_germ_start_where[0])
                                            + seq_before_j_num
                                            + j_line_pos[0]
                                            + 1
                                        ):
                                            j_cig_num = (
                                                int(
                                                    j_germ_finish_where[
                                                        (len(j_germ_finish_where) - 1)
                                                    ]
                                                )
                                                - int(j_germ_start_where[0])
                                                - cur_pos_is
                                                + seq_before_j_num
                                                + j_line_pos[0]
                                            )
                                            j_cig = j_cig + str(j_cig_num) + "="
                                        j_cigar = j_cig

                                    result["v_cigar"] = v_cigar
                                    result["d_cigar"] = d_cigar
                                    result["j_cigar"] = j_cigar
                                    result["V_mis_ratio"] = V_mis_ratio
                                    result["D_mis_ratio"] = D_mis_ratio
                                    result["J_mis_ratio"] = J_mis_ratio
                                    result["V_germline_mismatch"] = V_germline_mismatch
                                    result["D_germline_mismatch"] = D_germline_mismatch
                                    result["J_germline_mismatch"] = J_germline_mismatch
                                    result["v_germ_start_where"] = v_germ_start_where
                                    result["v_germ_finish_where"] = v_germ_finish_where
                                    result["d_germ_start_where"] = d_germ_start_where
                                    result["d_germ_finish_where"] = d_germ_finish_where
                                    result["j_germ_start_where"] = j_germ_start_where
                                    result["j_germ_finish_where"] = j_germ_finish_where
                                    result["v_cigar_mismatches"] = v_cigar_mismatches
                                    result["d_cigar_mismatches"] = d_cigar_mismatches
                                    result["j_cigar_mismatches"] = j_cigar_mismatches

                                elif (
                                    line[: len("Alignment summary")]
                                    == "Alignment summary"
                                ):
                                    totalDistanceFromVGene = 0
                                    totalAlignedLength = 0

                                    while True:
                                        next_line_sp = igb_handle.readline().split("\t")
                                        region = (
                                            next_line_sp[0].replace("-IMGT", "").lower()
                                        )
                                        if region == "total":
                                            break
                                        if region == "FR1":
                                            fr1_start = int(next_line_sp[1]) - 1
                                            fr1_end = int(next_line_sp[2])
                                            if sHelper.RepresentsInt(next_line_sp[6]):
                                                fr1_gap = int(next_line_sp[6])
                                            else:
                                                fr1_gap = 0
                                            while (
                                                fr1_end - fr1_start + fr1_gap
                                            ) % 3 != 0:
                                                fr1_start += 1
                                            result["fr1_from"] = fr1_start
                                            result["fr1_to"] = fr1_end
                                            result["fr1_length"] = int(next_line_sp[3])
                                            result["fr1_matches"] = int(next_line_sp[4])
                                            result["fr1_mismatches"] = int(
                                                next_line_sp[5]
                                            )
                                            result["fr1_gaps"] = int(next_line_sp[6])
                                            totalAlignedLength += result["fr1_length"]
                                            totalDistanceFromVGene += (
                                                result["fr1_mismatches"]
                                                + result["fr1_gaps"]
                                            )
                                        else:
                                            try:
                                                result[region + "_from"] = (
                                                    int(next_line_sp[1]) - 1
                                                )
                                                result[region + "_to"] = int(
                                                    next_line_sp[2]
                                                )
                                                result[region + "_length"] = int(
                                                    next_line_sp[3]
                                                )
                                                result[region + "_matches"] = int(
                                                    next_line_sp[4]
                                                )
                                                result[region + "_mismatches"] = int(
                                                    next_line_sp[5]
                                                )
                                                result[region + "_gaps"] = int(
                                                    next_line_sp[6]
                                                )
                                                totalAlignedLength += result[
                                                    region + "_length"
                                                ]
                                                totalDistanceFromVGene += (
                                                    result[region + "_mismatches"]
                                                    + result[region + "_gaps"]
                                                )
                                            except ValueError:
                                                continue

                                    result["v_alignment_length"] = totalAlignedLength
                                    result[
                                        "v_alignment_mutation"
                                    ] = totalDistanceFromVGene

                                elif (
                                    line[: len("V-(D)-J rearrangement summary")]
                                    == "V-(D)-J rearrangement summary"
                                ):
                                    next_line_sp = igb_handle.readline().split("\t")
                                    if chain_type in [
                                        ChainType.HUMAN_HEAVY,
                                        ChainType.RABBIT_HEAVY,
                                        ChainType.MOUSE_C57BL6_HEAVY,
                                        ChainType.MOUSE_BALBC_HEAVY,
                                        ChainType.MOUSE_HEAVY,
                                        ChainType.HUMAN_BETA,
                                        ChainType.HUMAN_DELTA,
                                        ChainType.MOUSE_BETA,
                                        ChainType.MOUSE_DELTA,
                                    ]:
                                        result["v_call"] = next_line_sp[0]
                                        result["d_call"] = next_line_sp[1]
                                        result["j_call"] = next_line_sp[2]
                                        result["vj_frame"] = next_line_sp[5].split("-")[
                                            0
                                        ]
                                        result["productive"] = next_line_sp[6]

                                    elif chain_type in [
                                        ChainType.HUMAN_LIGHT,
                                        ChainType.RABBIT_KAPPA,
                                        ChainType.MOUSE_C57BL6_LIGHT,
                                        ChainType.MOUSE_BALBC_LIGHT,
                                        ChainType.MOUSE_LIGHT,
                                        ChainType.HUMAN_ALPHA,
                                        ChainType.HUMAN_GAMMA,
                                        ChainType.MOUSE_ALPHA,
                                        ChainType.MOUSE_GAMMA,
                                    ]:
                                        result["v_call"] = next_line_sp[0]
                                        result["j_call"] = next_line_sp[1]
                                        result["vj_frame"] = next_line_sp[4].split("-")[
                                            0
                                        ]
                                        result["productive"] = next_line_sp[5]

                                # End of one alignment result
                                elif (
                                    line[: len("Effective search space used:")]
                                    == "Effective search space used:"
                                ):
                                    igblast_parsed_dict[query_id] = result
                                    # igblast_parsed_list.append(result)
                                    break

            elif domain_system == "kabat":
                print("Not supported domain_system...")
                return

        result_list = list(igblast_parsed_dict.values())
        return result_list

    def numbering_CDR3(
        self,
        chain_type: ChainType,
        cdr3_list: List[str],
        alphabet_with_number: bool = True,
    ) -> Tuple[TableHeader, TableData]:

        # check sequence type
        is_nucl = is_nucleotide_or_protein(cdr3_list[0])
        if is_nucl == None:
            raise RuntimeError("parameter cdr3_list unknown")

        if chain_type in (ChainType.HUMAN_LIGHT, ChainType.CHICKEN_LIGHT):
            CDR3_numbering = [
                "24",
                "25",
                "26",
                "27",
                "28",
                "29",
                "30",
                "A",
                "31",
                "32",
                "33",
                "34",
            ]
            len_end = 4
            cdr3_key = "LCDR3"
        elif chain_type in (ChainType.HUMAN_HEAVY, ChainType.CHICKEN_HEAVY):
            CDR3_numbering = ["95", "96", "97", "98", "99", "100", "A", "101", "102"]
            len_end = 2
            cdr3_key = "HCDR3"
        else:
            raise RuntimeError("parameter chain_type error")
        i_A = CDR3_numbering.index("A")
        if alphabet_with_number:
            last_num = CDR3_numbering[i_A - 1]
        else:
            last_num = ""
        CDR3_numbering[i_A] = last_num + CDR3_numbering[i_A]

        if is_nucl:
            max_len = int(max([len(s) / 3 for s in cdr3_list]))
        else:
            max_len = int(max([len(s) for s in cdr3_list]))
        init_len = len(CDR3_numbering)
        if max_len > init_len:
            CDR3_numbering = (
                CDR3_numbering[: i_A + 1]
                + [last_num + chr(ord("A") + i + 1) for i in range(max_len - init_len)]
                + CDR3_numbering[i_A + 1 :]
            )

        return_data = []
        for cdr3_str in cdr3_list:
            if cdr3_str == "-" or cdr3_str == "X" or cdr3_str == "":
                d = {cdr3_key: "-", "query": cdr3_str}
                for n in CDR3_numbering:
                    d[n] = ""
                return_data.append(d)
                continue

            if is_nucl:
                cdr3_str_aa = self.translate(cdr3_str)
            else:
                cdr3_str_aa = cdr3_str
            cdr3_len = len(cdr3_str_aa)
            d = {cdr3_key: cdr3_str_aa, "query": cdr3_str_aa}
            for n in CDR3_numbering:
                d[n] = ""

            if cdr3_len == 1:
                d[CDR3_numbering[0]] = cdr3_str_aa[0]
            else:
                d[CDR3_numbering[-1]] = cdr3_str_aa[-1]
                len_end_this = min(len_end, cdr3_len - 1)
                for i in range(-1, -1 * len_end_this - 1, -1):
                    d[CDR3_numbering[i]] = cdr3_str_aa[i]
                for i in range(cdr3_len - len_end_this):
                    d[CDR3_numbering[i]] = cdr3_str_aa[i]
            return_data.append(d)

        return CDR3_numbering, return_data

    def numbering_all(
        self,
        chain_type: ChainType,
        data: TableData,
        format_empty_residue="{bg_color:gray}",
    ) -> Tuple[TableHeader, List[TableHeader], List[TableData]]:

        if chain_type in (ChainType.HUMAN_LIGHT, ChainType.CHICKEN_LIGHT):
            numbering = {
                "FR1": [
                    "L0",
                    "L1",
                    "L2",
                    "L3",
                    "L4",
                    "L5",
                    "L6",
                    "L7",
                    "L8",
                    "L9",
                    "L10",
                    "L11",
                    "L12",
                    "L13",
                    "L14",
                    "L15",
                    "L16",
                    "L17",
                    "L18",
                    "L19",
                    "L20",
                    "L21",
                    "L22",
                    "L23",
                ],
                "CDR1": [
                    "L24",
                    "L25",
                    "L26",
                    "L27",
                    "L27A",
                    "L28",
                    "L29",
                    "L30",
                    "L31",
                    "L32",
                    "L33",
                    "L34",
                ],
                "FR2": [
                    "L35",
                    "L36",
                    "L37",
                    "L38",
                    "L39",
                    "L40",
                    "L41",
                    "L42",
                    "L43",
                    "L44",
                    "L45",
                    "L46",
                    "L47",
                    "L48",
                    "L49",
                ],
                "CDR2": ["L50", "L51", "L52", "L53", "L54", "L55", "L56"],
                "FR3": [
                    "L57",
                    "L58",
                    "L59",
                    "L60",
                    "L61",
                    "L62",
                    "L63",
                    "L64",
                    "L65",
                    "L66",
                    "L67",
                    "L68",
                    "L69",
                    "L70",
                    "L71",
                    "L72",
                    "L73",
                    "L74",
                    "L75",
                    "L76",
                    "L77",
                    "L78",
                    "L79",
                    "L80",
                    "L81",
                    "L82",
                    "L83",
                    "L84",
                    "L85",
                    "L86",
                    "L87",
                    "L88",
                ],
                "CDR3": [
                    "L89",
                    "L90",
                    "L91",
                    "L92",
                    "L93",
                    "L94",
                    "L95",
                    "L95A",
                    "L96",
                    "L97",
                ],
                "FR4": [
                    "L98",
                    "L99",
                    "L100",
                    "L101",
                    "L102",
                    "L103",
                    "L104",
                    "L105",
                    "L106",
                    "L106A",
                    "L107",
                    "L108",
                    "L109",
                ],
            }
        elif chain_type in (ChainType.HUMAN_HEAVY, ChainType.CHICKEN_HEAVY):
            numbering = {
                "FR1": [
                    "H0",
                    "H1",
                    "H2",
                    "H3",
                    "H4",
                    "H5",
                    "H6",
                    "H7",
                    "H8",
                    "H9",
                    "H10",
                    "H11",
                    "H12",
                    "H13",
                    "H14",
                    "H15",
                    "H16",
                    "H17",
                    "H18",
                    "H19",
                    "H20",
                    "H21",
                    "H22",
                    "H23",
                    "H24",
                    "H25",
                    "H26",
                    "H27",
                    "H28",
                    "H29",
                    "H30",
                ],
                "CDR1": ["H31", "H32", "H33", "H34", "H35", "H35A"],
                "FR2": [
                    "H36",
                    "H37",
                    "H38",
                    "H39",
                    "H40",
                    "H41",
                    "H42",
                    "H43",
                    "H44",
                    "H45",
                    "H46",
                    "H47",
                    "H48",
                    "H49",
                ],
                "CDR2": [
                    "H50",
                    "H51",
                    "H52",
                    "H52A",
                    "H53",
                    "H54",
                    "H55",
                    "H56",
                    "H57",
                    "H58",
                    "H59",
                    "H60",
                    "H61",
                    "H62",
                    "H63",
                    "H64",
                    "H65",
                ],
                "FR3": [
                    "H66",
                    "H67",
                    "H68",
                    "H69",
                    "H70",
                    "H71",
                    "H72",
                    "H73",
                    "H74",
                    "H75",
                    "H76",
                    "H77",
                    "H78",
                    "H79",
                    "H80",
                    "H81",
                    "H82",
                    "H82A",
                    "H83",
                    "H84",
                    "H85",
                    "H86",
                    "H87",
                    "H88",
                    "H89",
                    "H90",
                    "H91",
                    "H92",
                    "H93",
                    "H94",
                ],
                "CDR3": [
                    "H95",
                    "H96",
                    "H97",
                    "H98",
                    "H99",
                    "H100",
                    "H100A",
                    "H101",
                    "H102",
                ],
                "FR4": [
                    "H103",
                    "H104",
                    "H105",
                    "H106",
                    "H107",
                    "H108",
                    "H109",
                    "H110",
                    "H111",
                    "H112",
                    "H113",
                ],
            }
        else:
            raise RuntimeError("parameter chain_type error")

        regions = ["FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"]
        list_numbering_data = []
        list_numbering_header = [numbering[r] for r in regions]
        for i_r, (r, numbernig_header) in enumerate(
            zip(regions, list_numbering_header)
        ):
            numbering_data = []

            r_max_length = 0
            for d in data:
                d_num = {}
                if r not in d:
                    numbering_data.append(d_num)
                    continue

                is_nucl = is_nucleotide_or_protein(d[r])
                if is_nucl == None:
                    raise RuntimeError("unknown sequence type")
                if is_nucl:
                    r_seq_aa = self.translate(d[r])
                else:
                    r_seq_aa = d[r]
                r_max_length = max(r_max_length, len(r_seq_aa))

                last_index_alpha = -1
                for i in range(len(numbernig_header)):
                    if numbernig_header[-i - 1][-1].isalpha():
                        last_index_alpha = len(numbernig_header) - i - 1
                        break
                # if length is variable
                if r_max_length > len(numbernig_header):
                    if last_index_alpha != -1:
                        deficient_length = r_max_length - len(numbernig_header)
                        last_alpha = numbernig_header[last_index_alpha][-1]
                        numbernig_header = (
                            numbernig_header[: last_index_alpha + 1]
                            + [
                                numbernig_header[last_index_alpha][:-1]
                                + chr(ord(last_alpha) + j + 1)
                                for j in range(deficient_length)
                            ]
                            + numbernig_header[last_index_alpha + 1 :]
                        )
                        last_index_alpha += deficient_length
                        list_numbering_header[i_r] = numbernig_header

                if r == "FR1":
                    stick_left = 0
                elif r == "FR4":
                    stick_left = len(r_seq_aa)
                else:
                    # if do not have alphabet part
                    if last_index_alpha == -1:
                        stick_left = len(r_seq_aa)
                    else:
                        stick_left = len(r_seq_aa) - (
                            len(numbernig_header) - last_index_alpha - 1
                        )

                    if stick_left < 0:
                        stick_left = 0

                for i, n in enumerate(numbernig_header):
                    if i == stick_left:
                        break
                    d_num[n] = r_seq_aa[i]
                for j, n in enumerate(numbernig_header[::-1]):
                    if len(r_seq_aa) - stick_left == j:
                        break
                    d_num[n] = r_seq_aa[-j - 1]

                numbering_data.append(d_num)

            list_numbering_data.append(numbering_data)

        # fill the blank
        for r, numbernig_header, numbering_data in zip(
            regions, list_numbering_header, list_numbering_data
        ):
            for d_num in numbering_data:
                for n in numbernig_header:
                    if n not in d_num:
                        d_num[n] = format_empty_residue + "-"

        return regions, list_numbering_header, list_numbering_data

    def add_ptm_count_and_color(
        self, aa_seq_list=List[str], prefix_residue: str = "R", **kwargs
    ) -> Tuple[TableHeader, TableHeader, TableData]:
        """
        additional_ptm = [('D', '(D)'), ('E', '(E)'), ('D+E', '([DE])'), ('R', '(R)'), ('K', '(K)'), ('R+K', '([RK])')]

        :param header:
        :param data:
        :param col_aa_seq:
        :param col_residue: Typically, ('24', '34') for VL and ('95', '102') for VH.
                            If this value is None, just add count info and do not mark color.
        :return:
        """

        ret_data = []
        if "aa_residues" in kwargs and "resi_header" in kwargs:
            ret_data = kwargs["aa_residues"]
            resi_header = kwargs["resi_header"]
        else:
            lengths = [len(seq) for seq in aa_seq_list]
            max_length = int(max(lengths))
            for seq in aa_seq_list:
                d = {}
                d["query"] = seq
                for n in range(max_length):
                    if n < len(seq):
                        d[prefix_residue + str(n + 1)] = seq[n]
                    else:
                        d[prefix_residue + str(n + 1)] = ""
                ret_data.append(d)
            resi_header = [prefix_residue + str(n + 1) for n in range(max_length)]

        def mark_ptm(_d, _col_letters, _ptm_name, mark_str, match_list):
            _d[_ptm_name] = len(match_list)
            if len(_col_letters) > 0:
                for m in match_list:
                    for i in range(*m.span()):
                        _d[_col_letters[i]] = mark_str + _d[_col_letters[i]]

        default_ptms = [
            "Deamidation",
            "Aspartate_Isomerization",
            "N Glycosylation",
            "Cleavage",
            "Oxidation",
            "Free Cysteine",
            "Sulfation",
            "Methylation",
        ]

        finder = PTMFinder()
        for d in ret_data:
            seq = d["query"]
            col_letters = []
            for h in resi_header:
                try:
                    if d[h] != "" and detach_format(d[h]) == seq[len(col_letters)]:
                        col_letters.append(h)
                except Exception as ex:
                    print(ex)

            mark_ptm(
                d,
                col_letters,
                default_ptms[0],
                "{bg_color:blue}",
                finder.find_diamidation(seq),
            )
            mark_ptm(
                d,
                col_letters,
                default_ptms[1],
                "{bg_color:cyan}",
                finder.find_aspartate_isomerization(seq),
            )
            mark_ptm(
                d,
                col_letters,
                default_ptms[2],
                "{bg_color:gray}",
                finder.find_n_glycosylation(seq),
            )
            mark_ptm(
                d,
                col_letters,
                default_ptms[3],
                "{bg_color:green}",
                finder.find_cleavage(seq),
            )
            mark_ptm(
                d,
                col_letters,
                default_ptms[4],
                "{bg_color:lime}",
                finder.find_oxidation(seq),
            )
            mark_ptm(
                d,
                col_letters,
                default_ptms[5],
                "{bg_color:magenta}",
                finder.find_free_cysteine(seq),
            )
            mark_ptm(
                d,
                col_letters,
                default_ptms[6],
                "{bg_color:orange}",
                finder.find_sulfation(seq),
            )
            mark_ptm(
                d,
                col_letters,
                default_ptms[7],
                "{bg_color:yellow}",
                finder.find_methylation(seq),
            )

        return default_ptms, resi_header, ret_data

    def get_ptm_count(self, aa_seq_list=List[str]) -> Tuple[TableHeader, TableData]:
        default_ptms = [
            "Deamidation",
            "Aspartate_Isomerization",
            "N Glycosylation",
            "Cleavage",
            "Oxidation",
            "Free Cysteine",
            "Sulfation",
            "Methylation",
        ]

        ret_data = []
        finder = PTMFinder()
        for seq in aa_seq_list:
            d = {
                "query": seq,
                default_ptms[0]: len(finder.find_diamidation(seq)),
                default_ptms[1]: len(finder.find_aspartate_isomerization(seq)),
                default_ptms[2]: len(finder.find_n_glycosylation(seq)),
                default_ptms[3]: len(finder.find_cleavage(seq)),
                default_ptms[4]: len(finder.find_oxidation(seq)),
                default_ptms[5]: len(finder.find_free_cysteine(seq)),
                default_ptms[6]: len(finder.find_sulfation(seq)),
                default_ptms[7]: len(finder.find_methylation(seq)),
            }
            ret_data.append(d)

        return default_ptms, ret_data

    def get_pi_value(self, aa_seq: str) -> float:
        putil = ProteinAnalysis(aa_seq)
        return putil.isoelectric_point()

    def get_gravy_value(self, aa_seq: str) -> float:
        putil = ProteinAnalysis(aa_seq)
        return putil.gravy()

    def parse_IMPre(
        self,
        sequence_list: List[str],
        readcount_list: List[int],
        seq_type: str = "Ig",
        chain: str = "heavy",
        **kwargs
    ) -> TableData:
        """
        :param sequence_list:
        :param readcount_list:
        :param seq_type:
        :param chain:
        :return:
        """

        # check sequence type
        is_nucl = is_nucleotide_or_protein(sequence_list[0])
        if is_nucl == None or is_nucl == False:
            raise RuntimeError(
                "parameter sequence_list is not a nucleotide sequence list"
            )

        if "vn" in kwargs:
            vn = kwargs["vn"]
        else:
            vn = 300
        if "v_seed" in kwargs:
            v_seed = kwargs["v_seed"]
        else:
            v_seed = 200
        if "iter_num" in kwargs:
            iter_num = kwargs["iter_num"]
        else:
            iter_num = ""

        # csv to fa for run_IMPre().
        logging.info("---- Converting .csv to .fa started")
        tp_query = tempfile.NamedTemporaryFile("wt", suffix=".fa", delete=False)
        sequence_list_revised = np.repeat(np.array(sequence_list), readcount_list)
        for i, seq in enumerate(sequence_list_revised):
            tp_query.write(">seq_%d\n" % i)
            tp_query.write(seq + "\n")
        tp_query.close()
        logging.info("---- Converting .csv to .fa finished")

        logging.info("-- Running IMPre started")
        run_IMPre(
            query=tp_query.name,
            seq_type=seq_type,
            chain=chain,
            vn=vn,
            v_seed=v_seed,
            iter_num=iter_num,
        )
        logging.info("-- Running IMPre finished")

        result_dir = "/Tools/IMPre-master/Out_temp/%s_vn=%d_vseed=%d_iter=%d" % (
            chain,
            vn,
            v_seed,
            iter_num,
        )
        v_infer = open(
            os.path.join(result_dir, "%s_V_Germline.final.fa" % seq_type), "r"
        )
        j_infer = open(
            os.path.join(result_dir, "%s_J_Germline.final.fa" % seq_type), "r"
        )
        lines = v_infer.readlines()
        lines += j_infer.readlines()
        v_infer.close()
        j_infer.close()

        result_data = []
        for i in range(0, len(lines), 2):
            tmp_dict = {}
            tmp_dict["germline_ID"] = lines[i][1:-1]
            tmp_dict["sequence"] = lines[i + 1][:-1]
            result_data.append(tmp_dict)
        return result_data

    def blast_run(
        self,
        query_list: list(),
        db_list: list(),
        m_type: str,
        blast_dir: str = "/home/IP-team/data/blast_temp/",
        **kwargs
    ):
        if "file_blast" in kwargs:
            file_blast = kwargs["file_blast"]

        if "num_threads" in kwargs:
            num_threads = kwargs["num_threads"]
        else:
            num_threads = 10

        if "evalue" in kwargs:
            evalue = kwargs["evalue"]
        else:
            evalue = 0.1**5

        if "gapopen" in kwargs:
            gapopen = kwargs["gapopen"]
        else:
            gapopen = 1

        if "gapextend" in kwargs:
            gapextend = kwargs["gapextend"]
        else:
            gapextend = 2

        if "word_size" in kwargs:
            word_size = kwargs["word_size"]
        else:
            word_size = 10

        if "num_alignments" in kwargs:
            num_alignments = kwargs["num_alignments"]
        else:
            num_alignments = 1

        if "result_file_name" in kwargs:
            result_file_name = kwargs["result_file_name"]
        else:
            result_file_name = "blasted.txt"

        query_dir = os.path.join(blast_dir, "query")
        db_dir = os.path.join(blast_dir, "db")
        result_dir = os.path.join(blast_dir, "result")

        if not os.path.exists(query_dir):
            os.makedirs(query_dir)

        if not os.path.exists(db_dir):
            os.makedirs(db_dir)

        if not os.path.exists(result_dir):
            os.makedirs(result_dir)

        result_file = os.path.join(result_dir, result_file_name)

        with open(os.path.join(query_dir, "query.fasta"), "w") as fasta_writer:
            for i, e_seq in enumerate(query_list):
                fasta_writer.write(">" + str(i) + "\n")
                fasta_writer.write(e_seq + "\n")

        with open(os.path.join(db_dir, "db.fasta"), "w") as fasta_writer:
            for i, e_seq in enumerate(db_list):
                fasta_writer.write(">" + str(i) + "\n")
                fasta_writer.write(e_seq + "\n")

        if hostname == "JunhoLab":
            exe_path = "/Tools/ncbi-blast-2.7.1+/bin"
        elif hostname == "binel229":
            exe_path = "/home/team/IP-team/Tools/ncbi-blast-2.11.0+/bin"
        if m_type.upper() == "AA":
            # construct db
            format_cmd = "%s -in %s -dbtype prot -input_type fasta -out %s" % (
                os.path.join(exe_path, "makeblastdb"),
                os.path.join(db_dir, "db.fasta"),
                os.path.join(db_dir, "db"),
            )
            os.system(format_cmd)

            cline = NcbiblastpCommandline(
                num_threads=num_threads,
                query=os.path.join(query_dir, "query.fasta"),
                db=os.path.join(db_dir, "db"),
                out=result_file,
                num_descriptions=num_alignments,
                num_alignments=num_alignments,
                evalue=evalue,
            )

            format_exe = os.path.join(exe_path, "blastp")

            os.system(format_exe + str(cline)[len("blastp") :])

            print("Blastp completed.\n")

            handle = open(result_file, "r")
            blast_parser = NCBIStandalone.BlastParser()
            iterator = NCBIStandalone.Iterator(handle, blast_parser)

            resultList = []

            for i, each in enumerate(iterator):
                numAlignments = 0
                alignmentList = []

                if len(each.alignments) == 0:
                    tempSet = {
                        "num_alignments": numAlignments,
                        "alignments": alignmentList,
                    }
                    resultList.append(tempSet)
                else:
                    for alignment in each.alignments:
                        sbjctIndex = int(alignment.title[1:].strip())
                        hsps = alignment.hsps[0]
                        numAlignments += 1
                        queryStart = hsps.query_start
                        queryEnd = hsps.query_end
                        sbjctStart = hsps.sbjct_start
                        sbjctEnd = hsps.sbjct_end

                        matched = hsps.identities[0]
                        whole = hsps.identities[1]

                        numMismatches = hsps.identities[1] - hsps.identities[0]
                        mismatchList = []

                        matchPattern = ""

                        for i, eMatch in enumerate(hsps.match):
                            if eMatch == " ":
                                queryLoc = i + queryStart
                                sbjctLoc = i + sbjctStart
                                queryMismatch = hsps.query[i]
                                sbjctMismatch = hsps.sbjct[i]
                                mismatchList.append(
                                    {
                                        "query": queryMismatch,
                                        "sbjct": sbjctMismatch,
                                        "query_loc": queryLoc,
                                        "sbjct_loc": sbjctLoc,
                                    }
                                )
                                if queryMismatch == "-":
                                    matchPattern += "q"
                                elif sbjctMismatch == "-":
                                    matchPattern += "s"
                                else:
                                    matchPattern += " "
                            else:
                                matchPattern += "|"
                        alignmentList.append(
                            {
                                "matched": matched,
                                "whole": whole,
                                "sbjct_index": sbjctIndex,
                                "query_start": queryStart,
                                "query_end": queryEnd,
                                "sbjct_start": sbjctStart,
                                "sbjct_end": sbjctEnd,
                                "num_mismatches": numMismatches,
                                "mismatches": mismatchList,
                                "match_pattern": matchPattern,
                            }
                        )
                    tempSet = {
                        "num_alignments": numAlignments,
                        "alignments": alignmentList,
                    }
                    resultList.append(tempSet)
                return resultList

        if m_type.upper() == "NT":
            # construct db
            format_cmd = "%s -in %s -dbtype nucl -input_type fasta -out %s" % (
                os.path.join(exe_path, "makeblastdb"),
                os.path.join(db_dir, "db.fasta"),
                os.path.join(db_dir, "db"),
            )
            os.system(format_cmd)

            cline = NcbiblastnCommandline(
                num_threads=num_threads,
                query=os.path.join(query_dir, "query.fasta"),
                db=os.path.join(db_dir, "db"),
                evalue=evalue,
                out=result_file,
                gapopen=gapopen,
                gapextend=gapextend,
                word_size=word_size,
                num_descriptions=num_alignments,
                num_alignments=num_alignments,
            )

            format_exe = os.path.join("blastn")

            os.system(format_exe + str(cline)[len("blastn") :])

            print("Blastn completed.\n")

            handle = open(result_file, "r")
            blast_parser = NCBIStandalone.BlastParser()
            iterator = NCBIStandalone.Iterator(handle, blast_parser)

            resultList = []

            for i, each in enumerate(iterator):
                numAlignments = 0
                alignmentList = []

                if len(each.alignments) == 0:
                    tempSet = {
                        "num_alignments": numAlignments,
                        "alignments": alignmentList,
                    }
                    resultList.append(tempSet)
                else:
                    for alignment in each.alignments:
                        sbjctIndex = int(alignment.title[1:].strip())
                        hsps = alignment.hsps[0]
                        numAlignments += 1
                        queryStart = hsps.query_start
                        queryEnd = hsps.query_end
                        sbjctStart = hsps.sbjct_start
                        sbjctEnd = hsps.sbjct_end

                        matched = hsps.identities[0]
                        whole = hsps.identities[1]

                        numMismatches = hsps.identities[1] - hsps.identities[0]
                        mismatchList = []
                        matchPattern = ""
                        for i, eMatch in enumerate(hsps.match):
                            if eMatch == " ":
                                queryLoc = i + queryStart
                                sbjctLoc = i + sbjctStart
                                queryMismatch = hsps.query[i]
                                sbjctMismatch = hsps.sbjct[i]
                                mismatchList.append(
                                    {
                                        "query": queryMismatch,
                                        "sbjct": sbjctMismatch,
                                        "query_loc": queryLoc,
                                        "sbjct_loc": sbjctLoc,
                                    }
                                )
                                if queryMismatch == "-":
                                    matchPattern += "q"
                                elif sbjctMismatch == "-":
                                    matchPattern += "s"
                                else:
                                    matchPattern += " "
                            else:
                                matchPattern += "|"
                        alignmentList.append(
                            {
                                "matched": matched,
                                "whole": whole,
                                "sbjct_index": sbjctIndex,
                                "query_start": queryStart,
                                "query_end": queryEnd,
                                "sbjct_start": sbjctStart,
                                "sbjct_end": sbjctEnd,
                                "num_mismatches": numMismatches,
                                "mismatches": mismatchList,
                                "match_pattern": matchPattern,
                            }
                        )
                    tempSet = {
                        "num_alignments": numAlignments,
                        "alignments": alignmentList,
                    }
                    resultList.append(tempSet)
            return resultList


class PhyloWorker(Process):
    def __init__(self, in_queue, out_queue, run_type):
        Process.__init__(self)

        self.in_queue = in_queue
        self.out_queue = out_queue
        self.run_type = run_type

    def run(self):
        while True:
            # Get the work from the queue and expand the tuple
            target = self.in_queue.get()
            if target is None:
                self.in_queue.task_done()
                break
            if self.run_type == "extract_edges_phylotree":
                self.extract_edges_phylotree(**target)
            elif self.run_type == "build_tree":
                self.build_tree(**target)
            self.in_queue.task_done()

        logging.info("Finished run")

    def extract_edges_phylotree(
        self, col_group_id: str, group_name: str, msa_file: str
    ):
        tree_output = msa_file.replace("msa", "tree").replace(".fa", ".nex")
        is_exist = False
        if Path(tree_output).is_file() and os.stat(str(tree_output)).st_size > 0:
            tree = Phylo.read(tree_output, "nexus")
            is_exist = True
        else:
            align = AlignIO.read(msa_file, "fasta")
            calculator = DistanceCalculator("blastn")
            dm = calculator.get_distance(align)
            constructor = DistanceTreeConstructor()
            tree = constructor.upgma(dm)

        edges = []

        def get_clade_name(clade):
            return (
                group_name + "_" + clade.name if "Inner" in clade.name else clade.name
            )

        def get_edge_from_clade(clade, is_mother):
            if len(clade.clades) == 0:
                return
            edges.append(
                {
                    "from": get_clade_name(clade),
                    "to": get_clade_name(clade.clades[0]),
                    "length": clade.clades[0].branch_length,
                    "root": is_mother,
                    col_group_id: group_name,
                }
            )
            edges.append(
                {
                    "from": get_clade_name(clade),
                    "to": get_clade_name(clade.clades[1]),
                    "length": clade.clades[1].branch_length,
                    "root": 0,
                    col_group_id: group_name,
                }
            )
            get_edge_from_clade(clade.clades[0], 0)
            get_edge_from_clade(clade.clades[1], 0)

        get_edge_from_clade(tree.clade, 1)

        if not is_exist:
            Phylo.write(
                tree, msa_file.replace("msa", "tree").replace(".fa", ".nex"), "nexus"
            )

        self.out_queue.put([group_name, tree.clade.name, edges])

    def build_tree(self, aln_file: str, distance_model: str, **kwargs):
        align = AlignIO.read(aln_file, "fasta")
        calculator = DistanceCalculator(distance_model)
        dm = calculator.get_distance(align)
        constructor = DistanceTreeConstructor()
        tree = constructor.upgma(dm)
        self.out_queue.put({"tree": tree, **kwargs})

    @staticmethod
    def cut_tree(tree: Tree, branch_length_thres: float) -> List[Clade]:
        subtree_list = []

        def _dfs_cut_tree(_clade: Clade, _subtree_list: List[Clade]):
            new_children = []
            for i, _child in enumerate(_clade.clades):
                if _child.branch_length > branch_length_thres:
                    _child.branch_length = 0
                    _subtree_list.append(_child)
                else:
                    new_children.append(_child)
                _dfs_cut_tree(_child, _subtree_list)
            _clade.clades = new_children

        subtree_list.append(tree.clade)
        _dfs_cut_tree(tree.clade, subtree_list)

        return subtree_list

    @staticmethod
    def sort_sequence_by_phylo(sequence_batch, tree, out_file):
        map_reads = {}
        for rec in sequence_batch:
            map_reads[rec.id] = rec

        def dps_tree(clade, out_list):
            if len(clade.clades) == 0:
                out_list.append(clade.name)
            else:
                dps_tree(clade.clades[0], out_list)
                dps_tree(clade.clades[1], out_list)
            return

        sorted_id_list = []
        dps_tree(tree.clade, sorted_id_list)

        records_to_out = []
        for id in sorted_id_list:
            records_to_out.append(map_reads[id])

        SeqIO.write(records_to_out, out_file, "fasta")


def insert_into_pComb3x(full_vl_seq: str, full_vh_seq: str, linker_seq: str) -> str:
    vector_seq_file = r"C:\Users\hhjunny\Lib\Pcomb3X_Vector.fa"
    before_insert = "CCAGGCGGCC"
    after_insert = "ACTAGTGGCC"
    vector = SeqIO.read(vector_seq_file, "fasta")

    before_insert_full = str(vector.seq)[
        : str(vector.seq).index(before_insert) + len(before_insert)
    ]
    after_insert_full = str(vector.seq)[str(vector.seq).index(after_insert) :]
    return (
        before_insert_full + full_vl_seq + linker_seq + full_vh_seq + after_insert_full
    )
