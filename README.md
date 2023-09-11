

**Introduction to KFERQs**

KFERQs are short, linear motifs in protein sequences which are key elements in chaperone-mediated autophagy (CMA) and endosomal microautophagy (eMI). In both pathways the KFERQ motif is recognized by Hsp70 and its co-chaperones. After recognition the KFERQ-bearing protein is targeted for lysosomal degradation.
<img src="https://raw.githubusercontent.com/Marcin-Lubocki/KFERQminer/master/kferq_logo_ver_1.png" height="160" align="right" style="height:240px">


**The KFERQminer tool**

KFERQminer is an easy to use bioinformatics tool used for scanning user's protein sequences in order to find KFERQ-related motifs.
In order to support true positive hits we have developed the KFERQminer Database, which allows the exploration of already described KFERQ motifs. The collection is continuously updated and involves a thorough manual review of papers published in PubMed, searching for CMA-related keywords to extract and summarize the existing data. The collected set is subsequently used for sequence analysis and the development of new in silico KFERQ motif searching tool. 

You can find the collected motifs in **the KFERQminer Database** with annotated metadata here: https://github.com/Marcin-Lubocki/KFERQminer/blob/master/motif_db_expanded.csv and at the bottom of this notebook.

During the analysis your sequences will be checked for the presence of:
1. **hypothetical motifs** according to the common literature rules of KFERQ motif (aka canonical);
2. sequences **identical** to motifs collected from the literature;
3. sequences **similar** to motifs collected from the literature (it means they possess the same residue charge order).

In addition, we extract from the KFERQminer Database a set of high-quality motifs in order to build the first PSSM matrices, allowing us to assign mathematical scoring. To deal with the small dataset problem, we use: (i) Laplace smoothing to remove the zero probability problem, and (ii) an additional PSSM matrix based on a reduced alphabet by clustering amino acids based on their biochemical properties.

As the meta-analysis of the KFERQminer database showed that motifs requiring post-translational modifications are underrepresented, in order to support their further discovery we decided to use additional external software MusiteDeep (Wang et al., 2020) to predict whether a given motif can be phosphorylated or acetylated. In addition, using the knowledge that most short linear motifs are located in Intrinsically disordered regions, we calculate such probability for each motif using IUPred3 (Erdős et al., 2021). The whole pipeline is presented in Figure 2.

---


We wish to highlight some limitations:

- The records origin from manual PubMed review, so errors or missing motifs can occur. If you observe any error or missing record you can contact us via e-mail: marcin.lubocki@phdstud.ug.edu.pl. Thorough manual collection is very complex and time consumig process, though we will appreciate the help from other KFERQ researchers.
- As our goal was to collect the data published in peer-review journals, we did not question the quality of the results presented there. The metadata is therefore consistent with the authors' conclusions in the source publications.
- As a result, we included atypical motifs that do not meet the canonical rules, but which have been described in the literature.
- This is the beta version of the KFERQminer software. Please have in mind that you can approach a missed error. Please contact us if you find any.

---

**How to start**

Our tool has been made available as Jupyter Notebook. 
Also, to help non-technical users and ensure stability and computing power, we have made the entire software available in the form of a Jupyter Notebook via Google Colab (https://colab.research.google.com/drive/1QkC-qBJVku-iBEXDIQgkhjHUP-ZDZJOs?usp=sharing). Thanks to that the only requirement is to sign in to your google account.

In the Google Colab version, we have implemented comprehensive analysis by default using all available options. In this way, the user is only asked to load the input file, on which all analyzes will then be performed automatically. Please note that Google Colab has certain **limitations** regarding long analyzes, so for Big Data analyzes we suggest running the software in a local environment to avoid unnecessary errors. 
For non-technical people who encounter problems with running the tool in a local environment, **we will be happy to help by lending our computing resources**.
Moreover, for non-bioinformatics people, **we will be happy to help in performing more dedicated analyzes that will require dedicated development or tuning of the software**.

**To start your analysis on Google Colab just choose 'Executable environment' (at the very top of the page) and click 'Run all' or just press Ctrl+F9.**
You will be asked to upload your sequences in fasta file. Make sure they are in a proper format. The analysis will start automatically. After that just scroll down to see your results.

Empty tables mean your sequences did not have any KFERQ-related hits.

---

**External tools**

KFERQminer uses some external tools in its pipeline. The authors of KFERQminer are not the authors of these software. We therefore ask you to additionally cite the tools listed below.

**MusiteDeep**

https://www.musite.net/

If you find analysis related to the prediction of post-translational modification sites useful and would like to publish your data, **please cite the following papers for using MusiceDeep**:
Wang, D., et al. (2020) MusiteDeep: a deep-learning based webserver for protein post-translational modification site prediction and visualization, Nucleic Acids Research, Volume 48, Issue W1, 02 July 2020, Pages W140–W146.

Wang, D., et al. (2019) Capsule network for protein post-translational modification site prediction, Bioinformatics, 35(14), 2386-2394.

Wang, D., et al. (2017) MusiteDeep: a deep-learning framework for general and kinase-specific phosphorylation site prediction, Bioinformatics, 33(24), 3909-3916.


**IUPred3**

https://iupred3.elte.hu/

If you find analysis related to the intrinsically disordered region prediction useful and would like to publish your data, **please cite the following papers for using IUPred**:
Gábor Erdős, Mátyás Pajkos, Zsuzsanna Dosztányi

IUPred3: prediction of protein disorder enhanced with unambiguous experimental annotation and visualization of evolutionary conservation

Nucleic Acids Research 2021;49(W1):W297-W303.


Bálint Mészáros, Gábor Erdős, Zsuzsanna Dosztányi

IUPred2A: context-dependent prediction of protein disorder as a function of redox state and protein binding

Nucleic Acids Research 2018;46(W1):W329-W337.


Gábor Erdős, Zsuzsanna Dosztányi

Analyzing Protein Disorder with IUPred2A

Current Protocols in Bioinformatics 2020;70(1):e99


