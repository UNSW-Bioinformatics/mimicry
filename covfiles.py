
from pathlib import Path
import pandas as pd

app_dir = Path(__file__).parent
covlist = pd.read_csv(app_dir / "./coverage/totalcoverages.csv")
pcov = pd.read_csv(app_dir / "./coverage/indcov.csv")
world = pd.read_csv(app_dir / "./coverage/worldcov.csv")
regions = pd.read_csv(app_dir / "./coverage/regioncov.csv")
countries = pd.read_csv(app_dir / "./coverage/countrycov.csv")

def getPoptype(pop, regi, country, regions, countries):
    if regi != 'Select a region':
        regions = regions[regions['Region'].isin([regi])]
    if country != 'Select a country':
        countries = countries[countries['Country'].isin([country])]
    covlist = {'Worldwide' : (world, 'World'), 'by Region' : (regions, 'Region'), 'by Country' : (countries, 'Country')}
    return covlist[pop][0], covlist[pop][1]

def matchData(data, i, n):
    #'Cumulative Coverage (%)' : pcov[j][i]
    #print(pcov[str(n)][i])
    for j in range(1, n + 1):
        if str(covlist[f'Epitope {j}'][i]) == 'nan':
            break
        data.append({
            'No.' : j,
            'Allele' : covlist[f'HLA {j}'][i],
            'Best Epitope Match' : covlist[f'Epitope {j}'][i],
            'Net4 Rank (%)' : covlist[f'Net4 Rank {j}'][i],
            'Individual Coverage (%)' : covlist[f'Coverage % {j}'][i],
            'Cumulative Coverage (%)' : pcov[str(j)][i]
        })
    return pd.DataFrame(data)

def totalCov(pcov, cancer, area, ctype, retstr, num, i):
    s = '' if num == '1' else 's'
    if str(pcov[num][i]) != 'nan':
        return f'Predicted Coverage of {num} Epitope Combination{s} for {ctype} Cancer {retstr}: {str(pcov[num][i])}%.'
    j = int(num)
    while j > 0:
        if pcov['Cancer'][i] == cancer and pcov['Area'][i] == area and str(pcov[str(j)][i]) != 'nan':
            break
        j -= 1
    return f'{ctype} Cancer only has {j} Available Matches {retstr}, with a Predicted Coverage of {str(pcov[str(j)][i])}%.'

covstring = "This page shows the predicted coverage for an N-epitope cancer vaccine for different cancer types in different areas of the world, up to 30 combinations. Potential epitope matches are epitopes of the given cancer that also occur in the human microbiome. Using Net4 Rank percentages, the matching epitope with the strongest binding affinity to the given allele is shown. A lower Net4 Rank percentage is indicative of a stronger binding affinity of the epitope to the allele."
