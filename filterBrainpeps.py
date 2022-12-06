from filter import nmPepFilter
import pandas as pd
a='seq-1L'
base = pd.read_csv("data/new_dataset_sequence.csv")
naturals = []
largura_desejada=420

pd.set_option('display.width', largura_desejada)
pd.set_option('display.max_columns',10)

for x in base[a]:
    naturals.append(nmPepFilter(str(x)))

print(naturals)

result = pd.concat([base, pd.DataFrame(naturals)], axis=1)
result= result.rename(columns={0: 'NM'})
result.to_csv('results/dataset_sequence_NM.csv', index=False)