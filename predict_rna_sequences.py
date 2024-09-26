import sys
import os
import argparse
from Bio import SeqIO
import pandas as pd
from rna_analysis import *
from dna_analysis import *
import RNA_Inference

def read_fasta(file_path):
    return {record.id: str(record.seq) for record in SeqIO.parse(file_path, "fasta")}

def predict_seq(seq_id, full_seq, rna_inference, output_dir):
    target_columns = ['reactivity', 'deg_Mg_pH10', 'deg_pH10', 'deg_Mg_50C', 'deg_50C']
    
    predictions, attn_weight, input_features = rna_inference.predict(full_seq)
    
    rna_predictions_df = pd.DataFrame(columns=['position'] + target_columns)
    rna_predictions_df['position'] = np.arange(len(predictions))
    rna_predictions_df[target_columns] = predictions
    
    output_file = os.path.join(output_dir, f"{seq_id}_predictions.csv")
    rna_predictions_df.to_csv(output_file, index=False)
    print(f"Predictions for sequence {seq_id} saved to {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Predict RNA degradation with CDS + optional UTR.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file containing CDS sequences")
    parser.add_argument("-o", "--output", required=True, help="Output directory for prediction results")
    parser.add_argument("--utr-path", help="Input FASTA file containing UTR sequences (optional)")
    parser.add_argument("--model_weights", default="RNA_Inference/best_weights", help="Path to model weights (optional)")
    
    args = parser.parse_args()
    os.makedirs(args.output)

    # Load UTR (empty strings if not provided)
    utr_5, utr_3 = '', ''
    if args.utr_path:
        utr_sequences = read_fasta(args.utr_path)
        utr_5 = utr_sequences['u5:0']
        utr_3 = utr_sequences['u3:0']

    cds_sequences = read_fasta(args.input)

    # Load model
    rna_inference = RNA_Inference.RNA_Inference()
    rna_inference.load_models(args.model_weights)

    # Predict
    for cds_id, cds_seq in cds_sequences.items():
        full_seq = utr_5 + cds_seq + utr_3
        predict_seq(cds_id, full_seq, rna_inference, args.output)

if __name__ == "__main__":
    main()