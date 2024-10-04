Kaleb_Rotation ReadMe 

  IMPORTANT NOTE: Scripts will not work unless working directory is changed.

  Contains all information on data processing from raw counts to various plot analyses 

  The file "Protocol_Counts_Information" contains more detailed information about the accession number, protocol,     paper name, etc. 
  
  General Acronym Key Below 
  
    CN = Cortical Neurons 
    MG = Microglia 
    iPSC = Induced Pluripotent Stem Cell 
    DoD = Date of Differentiation 
    Proto = Protocols 
    
*Folders 

  Description of File Contents Listed Below 
  
-Folder 1: Combined Datasets 
  
  Contains all combined cleaned datasets. 
    
  Data was taken and processed manually from the "Raw Counts" folder. 
    
-Folder 2: Plots 

  Contains all plots obtained from running the various scripts in the "Scripts" folder 
  
  Plot Name Key Listed Below 
  
      Volc = Volcano Plots 
      HM = Heatmap 
      Bar = Bar Graphs 
      Enrich = Enrichment Analysis 
      
-Folder 3: Raw Counts 

  Contains All Raw Counts Obtained from GEO Database
  
  Tian Raw Counts isn't in the folder because the file was too large 
  
  GE Accession Numbers Listed Below 
  
      Ramadoss <- GSE272812
      Tian <- GSE124703
      Lisowski <- GSE233916
      Susco <- GSE144857
      Rummel <- GSE225813
      Hammonds <- GSE248998
      Breitmeyer <- GSE213232
      Penney <- GSE241858
      Guaqueta <- GSE221013
      
-Folder 4: Scripts

  Contains all scripts used to analyze the combined datasets

  Commands that are often edited for analysis are denoted by a "###" markdown intstead of a "#" markdown

  Script Function Key Below

    Analysis_Combination: Used to combine and clean manually processed raw counts into a combined datasets
    Analysis_DoD: Taken from Tian et al. Used to analyze differences in gene experession between dates of differentiation
    Analysis_Proto: All data is from CN samples. Used to analyze differences in gene experession between protocols
    Analysis_CN_iPSCs: Used to analye differences in gene experession between CN and iPSCs
    Analysis_MG_CN: Used to analyzed differnces in gene experession between MG and iPSCs
    Analysis_MG_CN_iPSC: Used to analyze differences in gene experession between MG, CN, and iPSCs

*Future Directions

  1. Date of differentiation and protocol analyses need to be performed for microglia cells.
  2. More cells types such as Adult and Fetal MG need to be added to analyses.
  3. Enirchment analyses need to be improved by selecting for smaller pathway groups.
  4. Enrichment tables need to be improved by adding a column depicting which pathway each gene belongs to.
