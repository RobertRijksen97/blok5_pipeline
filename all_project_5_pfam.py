# Datum: 24-10-2019
# Functie: Dit script is een pipline wat in staat is om een input FASTA bestand
# met MAFFT te multiple alignen, om hierna een HMM te bouwen en te zoeken met
# dit gebouwde HMM in de NR (gedownload) database
# Auteurs: Lars, Johannes en Robert

import os
import subprocess
from Bio import SearchIO, SeqIO


def main():
    """ Deze main functie roept alle functies aan en geeft doormiddel van een
    input de juist bestanden mee aan deze functies
    """
    try:
        naam_bestand_mafft_input = input("Van welk FASTA bestand wilt u een "
                                         "multiple alignment maken?: ")
        naam_output = input("Hoe moet het output bestand wilt u een multiple "
                            "alignment heten?: ")
        mafft(naam_bestand_mafft_input, naam_output)
        bestand_hmmbuild = input("Hoe wilt u het bestand noemen waarin de HMM "
                                 "gebouwd wordt: ")
        hmm_builder(bestand_hmmbuild, naam_bestand_mafft_input)
        database_nr = "/home/robert/Downloads/nr_database"
        accessions, descript = hmm_searcher(bestand_hmmbuild, database_nr)
        file_zetter(accessions, descript, database_nr)
    except FileNotFoundError:
        print("Bestand niet gevonden")
    except TypeError:
        print("Bestand niet gevonden --> check of het path naar de database"
              "klopt!")


def mafft(input_bestand, output_bestand):
    """Draai ClustalO op de linux command line

    Input:
    input_bestand - str, naam van fasta bestand om MSA van te maken
    output_bestand - str, naam van output van ClustalO

    Return:
    bestand wordt weggeschreven door ClustalO, deze kan in huidige path worden
    gevonden, tenzij anders aangegeven
    """
    try:
        # Check eerst of het outputbestand bestaat
        if os.path.isfile(output_bestand):
            pass
        # Mocht het bestand niet bestaan, dan voeten we het commando uit
        else:
            cmd = '"/usr/bin/mafft"  --retree 1 --reorder "{}" > "{}"' \
                .format(input_bestand, output_bestand)
            e = subprocess.check_call(cmd, shell=True)
            print(e)
        return output_bestand
    except subprocess.CalledProcessError:
        print("Fout in aanroep subprocess")


def hmm_builder(bestand_hmmbuild, bestand_hmminput):
    """ Deze functie bouwt de hmm

    :param bestand_hmmbuild: Hoe je het bestand wilt noemen waarin de hmm score
    staat
    :param bestand_hmminput: Waartegen je de hmm wilt bouwen
    """

    try:
        if os.path.isfile(bestand_hmmbuild):
            pass
        else:
            subprocess.run(["hm build", "{}.hmm".format(bestand_hmmbuild), "{}"
                           .format(bestand_hmminput)])
    except FileNotFoundError:
        print("Fout in subprocess of het bestand wordt niet gevonden")
    except subprocess.CalledProcessError:
        print("Fout in aanroep subprocess")


def hmm_searcher(bestand_hmm, bestand_zoeken):
    """ Deze functie zoekt mbv de hmm in het gegeven bestand

    :param bestand_hmm: het bestand waarin de hmm gebouwd is
    :param bestand_zoeken: Is het bestand waarin je de hmm wilt zoeken

    :return: Alle accessiecodes die een goede hit hebben
    """

    try:
        subprocess.run(['hmmsearch', '-o', 'output_hmm.txt', '{}.hmm'
                       .format(bestand_hmm), '{}'.format(bestand_zoeken)])
        accessions = []
        descript = []
        with open('output_hmm.txt') as intro:
            for qresult in SearchIO.parse(intro, 'hmmer3-text'):
                for hitresult in qresult:
                    id_hit = hitresult.id  # accessiecode met een hit
                    dest = hitresult.description  # naam / beschrijving met
                    # een hit
                    ev = hitresult.evalue
                    if ev < 0.01:  # alles kleiner dan een e-value van 0.1
                        # heeft een kans op homologe relatie
                        accessions.append(id_hit)
                        descript.append(dest)

        return accessions, descript

    except FileNotFoundError:
        print("Fout in subprocess of het bestand wordt niet gevonden")
    except subprocess.CalledProcessError:
        print("Fout in aanroep subprocess")
    except TypeError:
        print("Waarschijnlijke fout in bestandaanroep")


def file_zetter(acc, dest, bestand):
    """ Deze functie zoekt de sequenties bij de goed hmm search hits

    :param acc: Een lijst met accessiecodes met een goede hit
    :param dest: Een lijst met de beschrijvingen / naam met een goede hit
    :param bestand: Het bestand, waarin de hmmsearch is uitgevoerd, waarin
    gezocht kan worden naar de sequenties

    :return: Een bestand met alle resultaten netjes in een bestand geschreven
    welke in 1x in de database gezet kan worden
    """
    try:
        bestand2 = open("FINAL_RESULTS.txt", "a+")
        for seq_record in SeqIO.parse(bestand, "fasta"):
            # Omdat je de sequentie niet terug kan vragen bij de hmm search
            # was het noodzakelijk om deze functie te schrijven die de database
            # doorzoekt om de sequentie's te achterhalen met een juiste hit.
            for i in range(len(acc)):
                if acc[i] in seq_record.id:
                    bestand2.write(acc[i] + "\n")
                    bestand2.write(str(seq_record.seq) + "\n")
                    bestand2.write(dest[i] + "\n")
                    bestand2.write(str(1))
                    bestand2.write("\n")
                    bestand2.write(str(1))
                    bestand2.write("\n")

        bestand2.close()

    except FileNotFoundError:
        print("Bestand wordt niet gevonden")


main()
