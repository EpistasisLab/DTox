# DTox code structure

```
DTox code
|-- targettox.py: derive target binding profile of query compounds using the predictive MACCS fingerprints identified by TargetTox (feature selection pipeline)
|-- dtox.py: learn and evaluate DTox model
|   |-- dtox_data.py: contain data-formatting functions used for DTox model training.
|   |-- dtox_hierarchy.py: contain functions used to process sorted DTox hiearchy files and compute model statistics.
|   |-- dtox_nn.py: contain functions used to build basic neural network structure for DTox model.
|   |-- dtox_loss.py: contain the the loss function used in DTox model.
|   |-- early_stop.py: contain early stop function of DTox model.
|   |-- dtox_learning.py: contain deep learning functions used in the DTox model construction.
|-- dtox_interpret.py implement layer-wise relevance propagation to evaluate relevance of DTox paths.
|   |-- dtox_lrp.py: contain functions used for implementing LRP to evaluate relevance of DTox paths.
|   |-- dtox_visualize.R: use visNetwork package to visualize the flow of relevance along VNN paths between query compound, target protein, pathways, and the toxicity outcome.
```

