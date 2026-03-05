from pathlib import Path

import pandas as pd

app_dir = Path(__file__).parent
df = pd.read_csv(app_dir / "PanCancer.csv")
myStyle = pd.read_csv(app_dir / "styles.css")

# name formatting
pairs = {'Cervix' : 'Cervical', 'Ovary' : 'Ovarian', 'Pancreas' : 'Pancreatic', 'Uterus' : 'Uterine'}
