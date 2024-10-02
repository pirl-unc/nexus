import pandas as pd

if __name__ == '__main__':
    data = {
        'peptide': ['SIINFEKL',
                    'SIINFEKL'],
        'allele': ['H-2-Kb',
                   'H-2-Db']
    }
    pd.DataFrame(data).to_csv('../../../test/data/csv/sample900_mhcflurry2_input.csv', index=False)
