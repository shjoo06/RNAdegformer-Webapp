import sys
sys.path.append('draw_rna_pkg/')
import matplotlib.pyplot as plt
import os

from rna_analysis import *
from dna_analysis import *
import RNA_Inference

target_columns = [ 'reactivity', 'deg_Mg_pH10',
       'deg_pH10', 'deg_Mg_50C', 'deg_50C']
os.system('mkdir temp')
seq ="GGAAAAGCUCUAAUAACAGGAGACUAGGACUACGUAUUUCUAGGUAACUGGAAUAACCCAUACCAGCAGUUAGAGUUCGCUCUAACAAAAGAAACAACAACAACAAC"

# load model
rna_inference = RNA_Inference.RNA_Inference()
rna_inference.load_models('RNA_Inference/best_weights')

# predict, create df, save it to csv
predictions, attn_weight, input_features = rna_inference.predict(seq)
rna_predictions_df = pd.DataFrame(columns=['position']+target_columns)
rna_predictions_df['position']=np.arange(len(predictions))
rna_predictions_df[target_columns]=predictions
rna_predictions_df.to_csv('temp/rna_predictions.csv', index=False)
