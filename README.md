# digital_sorter
An R shiny application to identify cell surface markers by mimicking the FACS sorting process.

![image](https://github.com/user-attachments/assets/fcd09b23-b6a3-44ea-b0c5-6eb0acc3b925)
![image](https://github.com/user-attachments/assets/57f5e78c-f19c-48f0-8e2c-9deb209d8bfb)


# Developer guide
1. The _datasets.rds_ put in folder data/ should be a SplitObject of seurat.
   - Could use a separated R script **merge_seuratobject_addsplits_inspect.R** to create the split object (split by datasets) 
      and add **"split1", "split2", "split3", "split4"** based on the expression of master markers.
   - In the metadata, the annotation of cohort origins should be **"dataset_origin"**,
                       sample types should be **"disease"**.
   - Metadata should contain **"annotation.l1"** and **"annotation.l2"** from reference-based annotation (Azimuth)                
   
2. The _markerlist.rds_ in folder data/ should be a list including master markers of interested cancer types.

3. The _cell.surface.marker.rds_ in folder data/ should be a vector of cell surface markers.
   - Could be renew by function **get_cell_markers()** in **util.R**
