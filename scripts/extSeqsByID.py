import argparse
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment


def extraer_secuencias(msa_file, id_file, output_file):
    # Leer el archivo MSA
    alignment = AlignIO.read(msa_file, 'fasta')
    
    # Leer los IDs de las secuencias a extraer
    with open(id_file, 'r') as f:
        ids = set([line.strip() for line in f])
    
    # Filtrar el MSA para contener solo las secuencias con IDs en la lista
    new_alignment = MultipleSeqAlignment([record for record in alignment if record.id in ids])
    
    # Escribir el nuevo MSA al archivo de salida
    AlignIO.write(new_alignment, output_file, 'fasta')
    print(f"Nuevo archivo MSA creado con {len(new_alignment)} secuencias: {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Extraer secuencias de un archivo MSA según una lista de IDs y crear un nuevo archivo MSA.')
    
    parser.add_argument('-m', '--msa_file', type=str, required=True, help='Archivo de entrada MSA en formato FASTA.')
    parser.add_argument('-i', '--id_file', type=str, required=True, help='Archivo de texto con los IDs de las secuencias a extraer.')
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Archivo de salida MSA en formato FASTA.')
    
    args = parser.parse_args()
    
    extraer_secuencias(args.msa_file, args.id_file, args.output_file)

if __name__ == "__main__":
    main()

