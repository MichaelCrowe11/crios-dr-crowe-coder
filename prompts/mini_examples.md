# Few-Shot Mini Examples

## 1. Filtering & Ranking (Python / R)
```python
# Python (pip install rdkit-pypi)
from rdkit import Chem
from rdkit.Chem import Descriptors, QED

def filter_rank(smiles_list, mw_max=500, qed_min=0.55):
    out = []
    for s in smiles_list:
        m = Chem.MolFromSmiles(s)
        if not m:
            continue
        mw = Descriptors.MolWt(m)
        q = QED.qed(m)
        if mw <= mw_max and q >= qed_min:
            out.append((s, q))
    return sorted(out, key=lambda x: x[1], reverse=True)
```
```r
# R (install.packages('rcdk'); needs Java)
library(rcdk)
filter_rank <- function(smiles, mw_max=500, qed_min=0.55){
  res <- list()
  for (s in smiles){
    m <- parse.smiles(s)[[1]]
    if (is.null(m)) next
    # Approximate QED placeholder; real QED requires extra calc
    mw <- get.total.mass(m)
    q  <- 0.5
    if (mw <= mw_max && q >= qed_min) res[[length(res)+1]] <- list(s,q)
  }
  ord <- order(sapply(res, function(x) x[[2]]), decreasing=TRUE)
  res[ord]
}
```

## 2. Simple QSAR Baseline (Scaffold Split)
```python
# Python baseline (pip install rdkit-pypi scikit-learn)
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score
import random

def morgan_fp(smi):
    m = Chem.MolFromSmiles(smi)
    return AllChem.GetMorganFingerprintAsBitVect(m, 2, 2048) if m else None

# Placeholder dataset
smiles = ["CCO", "CCN", "c1ccccc1", "CC(=O)O"]
labels = [0,1,0,1]
X = [morgan_fp(s) for s in smiles]
X = [x for x in X if x is not None]
Y = labels
clf = RandomForestClassifier(n_estimators=64, random_state=42)
# Convert bitvectors to list features
import numpy as np
X_arr = np.array([list(x) for x in X])
clf.fit(X_arr, Y)
print("AUC", roc_auc_score(Y, clf.predict_proba(X_arr)[:,1]))
```
```r
# R baseline (install.packages(c('randomForest')) ; for real FP use rcdk)
library(randomForest)
smiles <- c('CCO','CCN','c1ccccc1','CC(=O)O')
labels <- c(0,1,0,1)
# Placeholder binary features (string length parity)
X <- data.frame(f1 = nchar(smiles) %% 2, f2 = nchar(smiles))
model <- randomForest(x=X, y=as.factor(labels), ntree=64)
print(model)
```

## Disclaimer
**Research planning only. No clinical advice or efficacy claims; verification and regulatory review required before any real-world use.**
