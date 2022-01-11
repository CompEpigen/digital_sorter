# digital_sorter
### 1. The **datasets.rds** put in folder data/ should be a SplitObject of seurat.
    - Could use a separated R script 20211208_merge_seuratobject_addsplits_inspect.R to create the split object (split by datasets) 
      and add "split1","split2","split3", "split4" based on the expression of master markers.
    - In the metadata, the annotation of cohort origins should be "dataset_origin",
                       sample types should be "disease".
    - Metadata should contain "annotation.l1" and "annotation.l2" from reference-based annotation (Azimuth)                
   
### 2. The **markerlist.rds** in folder data/ should be a list including master markers of interested cancer types.

### 3. The **cell.surface.marker.rds** in folder data/ should be a vector of cell surface markers.
    - Could be renew by function get_cell_markers() in util.R
