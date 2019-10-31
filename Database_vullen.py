import mysql.connector


def main():
    msa_seqs = open_msa_file()
    blast_seqs = open_protein_data_file()
    to_database(msa_seqs, blast_seqs)


def open_msa_file():
    """Reads through a MSA file and collects the headers and sequences

    Return:
    A dictionary containing the headers and sequences from a MSA file
    """
    file = open('final_results_mafft_output.txt', 'r')
    msa_seqs = {}
    header = ""
    seq = ""

    for line in file:
        if line.startswith(">"):
            msa_seqs[header] = seq
            seq = ""
            header = line.strip('\n')
        else:
            seq = seq + line.strip('\n')

    msa_seqs[header] = seq
    msa_seqs.pop("")
    print("klaar msa")
    return msa_seqs


def open_protein_data_file():
    """Reads through a file containing BLAST results and collects multiple parameters from this file

    Return:
    A dictionary containing multiple parameters from a BLAST results
    """
    file = open("FINAL_RESULTS22.txt", "r")
    blast_data = {}
    parameters = []
    accessiecode = ""

    for line in file:
        blast_data[accessiecode] = parameters
        parameters = []
        accessiecode = line.strip('\n')
        for _ in range(4):
            parameters.append(file.readline().strip('\n'))

    blast_data[accessiecode] = parameters
    blast_data.pop("")
    print("klaar sequentie data lezen")
    return blast_data


def to_database(msa_seqs, eiwit_seqs):
    """Puts results from an experiment in a database

    Input:
    msa_seqs:   A dictionary containing the headers and sequences from a MSA file
    blast_data: A dictionary containing multiple parameters from a BLAST results
    """
    print("begin vullen database")
    connection = mysql.connector.connect(host="hannl-hlo-bioinformatica-mysqlsrv.mysql.database.azure.com",
                                         user="sneoz@hannl-hlo-bioinformatica-mysqlsrv",
                                         db="sneoz",
                                         passwd="590830",
                                         port="3306")

    cursor = connection.cursor()

    cursor.execute("""insert into Eiwit_familie (
    Eiwit_fam_ID, Eiwit_familie) 
    values ({}, '{}')""".format(
        1, "HslJ-like superfamily"))

    cursor.execute("""insert into Onderzoeks_sequentie (
    Onderzoeks_sequentie, Naam_eiwit, Accessiecode, Eiwit_fam_ID, Onderzoek_sequentie_ID)
    values ('{}', '{}', '{}', {}, {})""".format(
        "MKKVAAFVALSLLMAGCVSNDKIAVTPEQLQHHRFVLESVNGKPVTSDKNPPEISFGEKMMISGSMCNRF"
        "SGEGKLSNGELTAKGLAMTRMMCANPQLNELDNTISEMLKEGAQVDLTANQLTLATAKQTLTYKLADLMN",
        "Heat shock protein HslJ", "P52644", 1, 1))

    for accessiecode, parameters in eiwit_seqs.items():
        cursor.execute("""insert into Eiwit_sequenties (
        Sequentie, Naam_eiwit, Accessiecode, Onderzoek_sequentie_ID, Generatie)
         values ('{}', '{}', '{}', {}, {})""".format(
            parameters[0], parameters[1], accessiecode, parameters[2], parameters[3]))

    cursor.execute("""SET FOREIGN_KEY_CHECKS = 0""")
    cursor.execute("""truncate table MSA_MAFFT""")

    for header, seq in msa_seqs.items():
        cursor.execute("""insert into MSA_MAFFT (
        Sequentie_aligned, Accessiecode) 
        values ('{}', '{}')""".format(
            seq, header))

    cursor.execute("""SET FOREIGN_KEY_CHECKS = 1""")

    connection.commit()
    print("helemaal klaar")


main()
