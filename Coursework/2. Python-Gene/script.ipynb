{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "ea5404fa-1d7e-4a8d-ad98-e8ec4925b5e0",
   "metadata": {},
   "outputs": [],
   "source": [
    "#Summary Statistics:\n",
    "!pip install pandas"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 33,
   "id": "5d1ca8bb-2b30-4c3b-ba91-affd3869abf5",
   "metadata": {},
   "outputs": [],
   "source": [
    "import pandas as pd\n",
    "df_KF= pd.read_csv(\"GeneExpression/Datasets/set_16/k_vs_F.deseq2.results.tsv\", sep = \"\\t\")\n",
    "df_KL= pd.read_csv(\"GeneExpression/Datasets/set_16/k_vs_L.deseq2.results.tsv\", sep = \"\\t\")"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "8b46e974-a393-4a8a-9cc6-117d336a8ee3",
   "metadata": {},
   "outputs": [],
   "source": [
    "df_KF.head()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "7c984e8f-b913-4bc6-96f5-7e0e06f4fa54",
   "metadata": {},
   "outputs": [],
   "source": [
    "df_KL.head()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 20,
   "id": "f4079986-37a3-438b-9fb7-6bd1b0e30ad0",
   "metadata": {},
   "outputs": [],
   "source": [
    "padj_cutoff = 0.05 # Adjusted p-value threshold (significance requires padj < this value)\n",
    "logFC_up = 1  # Threshold for upregulation: log2FoldChange > 1\n",
    "logFC_down = -1"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 21,
   "id": "2d0aa9f5-d03e-47ce-8b14-8e52d2a811d5",
   "metadata": {},
   "outputs": [],
   "source": [
    "KF_up = df_KF[(df_KF[\"padj\"] < padj_cutoff) & (df_KF[\"log2FoldChange\"] > logFC_up)] #Select significantly upregulated genes (padj < 0.05 AND log2FC > 1)\n",
    "KF_down = df_KF[(df_KF[\"padj\"] < padj_cutoff) & (df_KF[\"log2FoldChange\"] < logFC_down)] #Select significantly downregulated genes (padj < 0.05 AND log2FC < -1)\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 22,
   "id": "58818cd2-c9c4-437e-8aca-2ee64b803d59",
   "metadata": {},
   "outputs": [],
   "source": [
    "KL_up = df_KL[(df_KL[\"padj\"] < padj_cutoff) & (df_KL[\"log2FoldChange\"] > logFC_up)]\n",
    "KL_down = df_KL[(df_KL[\"padj\"] < padj_cutoff) & (df_KL[\"log2FoldChange\"] < logFC_down)]\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "7542da32-fbd4-4c97-8968-038c97d3bc9a",
   "metadata": {},
   "outputs": [],
   "source": [
    "summary = pd.DataFrame({\n",
    "    \"Comparison\": [\"K_vs_F\", \"K_vs_L\"],\n",
    "    \"Upregulated\": [len(KF_up), len(KL_up)],\n",
    "    \"Downregulated\": [len(KF_down), len(KL_down)],\n",
    "    \"Total_DEGs\": [len(KF_up)+len(KF_down), len(KL_up)+len(KL_down)]\n",
    "})\n",
    "\n",
    "summary\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "af5635a0-25b0-4b21-91c9-6ef436aa7b26",
   "metadata": {},
   "outputs": [],
   "source": [
    "#Summary of p-values and log fold changes across all genes for each comparison.\n",
    "def summary_stats(df):\n",
    "    print(\"\\nSummary of p-values:\")\n",
    "    print(\"Mean:\", df[\"pvalue\"].mean())\n",
    "    print(\"Median:\", df[\"pvalue\"].median())\n",
    "    print(\"Min:\", df[\"pvalue\"].min())\n",
    "    print(\"Max:\", df[\"pvalue\"].max())\n",
    "\n",
    "    print(\"\\nSummary of log2FoldChange:\")\n",
    "    print(\"Mean:\", df[\"log2FoldChange\"].mean())\n",
    "    print(\"Median:\", df[\"log2FoldChange\"].median())\n",
    "    print(\"Min:\", df[\"log2FoldChange\"].min())\n",
    "    print(\"Max:\", df[\"log2FoldChange\"].max())\n",
    "\n",
    "\n",
    "print(\"\\n=== K vs F Summary ===\")\n",
    "summary_stats(df_KF)\n",
    "\n",
    "print(\"\\n=== K vs L Summary ===\")\n",
    "summary_stats(df_KL)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "9f1c32e7-a88d-49e9-9153-8068e37806d6",
   "metadata": {},
   "outputs": [],
   "source": [
    "!pip install matplotlib "
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "580a2009-664e-4bc3-a76c-cda42c272c7e",
   "metadata": {},
   "outputs": [],
   "source": [
    "import seaborn as sns # Import seaborn for  scatter plotting\n",
    "import numpy as np # Import numpy for numerical operations\n",
    "import matplotlib.pyplot as plt # Import matplotlib for plotting functions\n",
    "\n",
    "def Volcano_plot(df, title, padj_cutoff=0.05, logFC_cutoff=1):\n",
    "    df = df.copy()\n",
    "    df[\"-log10padj\"] = -np.log10(df[\"padj\"])   # Calculation of Significance on the Y-axis\n",
    "    df[\"status\"] = \"Not Sig\" #gene classification\n",
    "    df.loc[(df[\"padj\"]<padj_cutoff) & (df[\"log2FoldChange\"]>logFC_cutoff),\"status\"] = \"Up\"\n",
    "    #Genes with adjusted p-value < 0.05 and log2FC > 1 were classified as Up-regulated,while those with log2FC < −1 were defined as Down-regulated.\n",
    "    df.loc[(df[\"padj\"]<padj_cutoff) & (df[\"log2FoldChange\"]<-logFC_cutoff),\"status\"] = \"Down\"\n",
    "    colors = {\"Up\":\"red\", \"Down\":\"blue\", \"Not Sig\":\"gray\"}\n",
    "    plt.figure(figsize=(6,5))\n",
    "    \n",
    "    for group in [\"Not Sig\",\"Up\",\"Down\"]:\n",
    "        sub = df[df[\"status\"]==group]\n",
    "        plt.scatter(sub[\"log2FoldChange\"],sub[\"-log10padj\"],\n",
    "                    color=colors[group],s=15,label=group)\n",
    "        \n",
    "    plt.axvline(logFC_cutoff, color=\"black\", linestyle=\"--\")\n",
    "    plt.axvline(-logFC_cutoff, color=\"black\", linestyle=\"--\")\n",
    "    plt.axhline(-np.log10(padj_cutoff), color=\"black\", linestyle=\"--\")\n",
    "\n",
    "    plt.xlabel(\"log2FoldChange\")\n",
    "    plt.ylabel(\"-log10(padj)\")\n",
    "    plt.title(title)\n",
    "    plt.legend(title=\"Differential Expression\")\n",
    "    plt.show()\n",
    "\n",
    "Volcano_plot(df_KF, \"Volcano_plot - KF\")\n",
    "Volcano_plot(df_KL, \"Volcano_plot - KL\")\n",
    "\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "df5522f5-c669-42ad-9197-08b011250715",
   "metadata": {},
   "outputs": [],
   "source": [
    "\n",
    "import pandas as pd\n",
    "import numpy as np\n",
    "import seaborn as sns\n",
    "import matplotlib.pyplot as plt\n",
    "\n",
    "def MA_plot(df, title, padj_cutoff=0.05, logFC_cutoff=1):\n",
    "    df = df.copy()\n",
    "    df[\"logBaseMean\"] = np.log10(df[\"baseMean\"] + 1)  # Convert mean expression to log10 scale for MA plot  \n",
    "    #gene classification\n",
    "    df[\"status\"] = \"Not Sig\"\n",
    "    df.loc[(df[\"padj\"]<padj_cutoff) & (df[\"log2FoldChange\"]>logFC_cutoff),\"status\"] = \"Up\"\n",
    "    df.loc[(df[\"padj\"]<padj_cutoff) & (df[\"log2FoldChange\"]<-logFC_cutoff),\"status\"] = \"Down\"\n",
    "    colors = {\"Up\":\"red\", \"Down\":\"blue\", \"Not Sig\":\"gray\"}\n",
    "    plt.figure(figsize=(6,5))\n",
    "    \n",
    "    for group in [\"Not Sig\",\"Up\",\"Down\"]:\n",
    "        sub = df[df[\"status\"]==group] #Filter dataframe to only include genes of this category\n",
    "        plt.scatter(sub[\"logBaseMean\"], sub[\"log2FoldChange\"],\n",
    "                    s=15, color=colors[group], label=group)\n",
    "        \n",
    "    plt.axhline(0, color=\"black\", linestyle=\"--\") #Horizontal reference line showing no fold-change boundary\n",
    "    plt.xlabel(\"log10(baseMean)\")\n",
    "    plt.ylabel(\"log2FoldChange\")\n",
    "    plt.title(title)\n",
    "    plt.legend(title=\"Expression Change\")\n",
    "    plt.show()\n",
    "\n",
    "\n",
    "MA_plot(df_KF, \"MA Plot - KF\")\n",
    "MA_plot(df_KL, \"MA Plot - KL\")"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "1e1cbb4d-d586-4b35-874a-fa79b3c81501",
   "metadata": {},
   "outputs": [],
   "source": [
    "#Histogram of p-values to assess the distribution of statistical significance.\n",
    "import pandas as pd\n",
    "import numpy as np\n",
    "import seaborn as sns\n",
    "import matplotlib.pyplot as plt\n",
    "\n",
    "def pval_hist(df, title):\n",
    "    plt.figure(figsize=(6,5))\n",
    "    plt.hist(df[\"pvalue\"], bins=50, color=\"skyblue\", edgecolor=\"black\")\n",
    "    plt.xlabel(\"p-value\")\n",
    "    plt.ylabel(\"Frequency\")\n",
    "    plt.title(title)\n",
    "    plt.show()\n",
    "\n",
    "pval_hist(df_KF, \"P-value Histogram - KF\")\n",
    "pval_hist(df_KL, \"P-value Histogram - KL\")\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "4007d9c2-7f9b-442c-813e-78e33eede8d2",
   "metadata": {},
   "outputs": [],
   "source": [
    "#Heatmap of the top differentially expressed genes to illustrate gene expression patterns across the conditions.\n",
    "\n",
    "import pandas as pd\n",
    "import seaborn as sns\n",
    "import matplotlib.pyplot as plt\n",
    "\n",
    "def heatmap_DEG_logFC(df, title, topN=30):\n",
    "   #  Select significantly differentially expressed genes\n",
    "    df_sig = df[df[\"padj\"] < 0.05].copy() ##  Select significantly differentially expressed genes\n",
    "    df_sig = df_sig.sort_values(\"log2FoldChange\", key=abs, ascending=False).head(topN) #sort genes by absolute fold change magnitude,keep Top N strongest DEGs\n",
    "    heatmap_data = df_sig[[\"log2FoldChange\"]] #Select only \"log2FoldChange\" column for heatmap\n",
    "    heatmap_data.index = df_sig[\"gene_id\"] #Set \"gene_id\" as row index to label genes in heatmap\n",
    "    plt.figure(figsize=(10,14))\n",
    "    sns.heatmap(heatmap_data, cmap=\"bwr\", center=0, annot=True, fmt=\".2f\") #Set the image size, draw a heatmap, with red indicating an increase and blue indicating a decrease, with 0 as the color center, and keep the precision to two decimal places.\n",
    "    plt.title(f\"{title} — Top {topN} DEGs\")\n",
    "    plt.xlabel(\"log2FoldChange (Expression Trend)\")\n",
    "    plt.ylabel(\"Genes\")\n",
    "    plt.savefig(f\"/Users/nicole/Desktop/{title}.png\", dpi=400)\n",
    "    plt.show()\n",
    "\n",
    "heatmap_DEG_logFC(df_KF,\"KF Heatmap\")\n",
    "heatmap_DEG_logFC(df_KL,\"KL Heatmap\")\n",
    "\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "ce3b7593-133f-4de7-ad38-93e76698fa8f",
   "metadata": {},
   "outputs": [],
   "source": [
    "#A table or list of significantly upregulated and downregulated genes with their corresponding fold changes, p-values, and adjusted p-values.\n",
    "def DEG_table(df, padj_cutoff=0.05, logFC_cutoff=1):\n",
    "    df = df.copy()\n",
    "    up = df[(df[\"padj\"] < padj_cutoff) & (df[\"log2FoldChange\"] > logFC_cutoff)]# Significantly upregulated genes\n",
    "    down = df[(df[\"padj\"] < padj_cutoff) & (df[\"log2FoldChange\"] < -logFC_cutoff)]# Significantly downregulated genes\n",
    "    print(f\" Significantly Upregulated genes: {len(up)}\")\n",
    "    print(f\" Significantly Downregulated genes: {len(down)}\\n\")\n",
    "    up_table = up[[\"gene_id\",\"log2FoldChange\", \"pvalue\", \"padj\"]]# Extracts columns of interest from significantly upregulated genes.\n",
    "    down_table = down[[\"gene_id\",\"log2FoldChange\", \"pvalue\", \"padj\"]]\n",
    "\n",
    "    return up_table, down_table\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 60,
   "id": "dab49092-09c5-4e99-a3b1-274cc8516c57",
   "metadata": {},
   "outputs": [],
   "source": [
    "up_KF, down_KF = DEG_table(df_KF)\n",
    "up_KL, down_KL = DEG_table(df_KL)\n",
    "\n",
    "up_KF.head(), down_KF.head()\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "7309690d-26e5-4ed0-81f1-6774c2128152",
   "metadata": {},
   "outputs": [],
   "source": []
  },
  {
   "cell_type": "code",
   "execution_count": 98,
   "id": "720a7210-3a71-4a43-ad6f-98a5d36ba6f1",
   "metadata": {},
   "outputs": [],
   "source": [
    "# Additional Analyses:\n",
    "Based solely on visual output, the experimental condition triggers non-random gene expression shifts, resulting in clear separation of regulated genes.\n",
    "Clustering patterns indicate that DEGs are functionally connected rather than independent responders.\n",
    "The system likely undergoes pathway-level regulation, where groups of genes are switched ON/OFF to adapt to the condition."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "c6c03935-807b-4d66-9d79-412d1b3dba9e",
   "metadata": {},
   "outputs": [],
   "source": [
    "\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "9b282c64-45e5-41ea-a305-0b2a5f329876",
   "metadata": {},
   "outputs": [],
   "source": []
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "f96573ce-8b11-4989-bda1-bc4ecfb46fb2",
   "metadata": {},
   "outputs": [],
   "source": []
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3 (ipykernel)",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.11.14"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}
