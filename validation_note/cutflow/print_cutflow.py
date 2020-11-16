import pandas as pd

df = pd.read_excel('cutflow.xlsx')
df = df.iloc[0:61,0:5]

#print(df.to_latex(index=True, escape=False, float_format = lambda x: '{:0.2f}'.format(x) if pd.notna(x) else '-'))

tab = df.to_latex(index=True, escape=False, float_format = lambda x: '{:0.2f}'.format(x) if pd.notna(x) else '-')
print tab
#print tab.replace('\\\\\n', '\\\\ \\hline\n')
