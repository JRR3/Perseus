# Import modules for scientific computing
# and sc-RNA-seq

import os
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sb
import matplotlib as mpl
import matplotlib.pyplot as plt


#The "frontend" is the user facing code, i.e., the plotting code,
#whereas the "backend" does all the hard work behind-the-scenes to
#make the figure. There are two types of backends: user interface
#backends (for use in PyQt/PySide, PyGObject, Tkinter, wxPython,
#or macOS/Cocoa); also referred to as "interactive backends")
#and hardcopy backends to make image files
#(PNG, SVG, PDF, PS; also referred to as 
#"non-interactive backends").

#Since we don't have a screen, we select the
#agg backend.
mpl.use('agg')


class Perseus:
    #This is a pipeline to do single-cell RNA sequencing.

    def __init__(self):

        print('Perseus is getting ready...')

        # What version of the libraries do we have?
        sc.logging.print_header()

        #Plotting settings
        sc.settings.set_figure_params(dpi=300,
                                      dpi_save=300,
                                      vector_friendly=True,
                                      facecolor='white')

        #Note that we do not have an interactive session.
        self.show_plot = False


        #Show all output, but make no suggestions.
        sc.settings.verbosity = 3


        mpl.rcParams['figure.dpi']=300
        #mpl.rcParams['figure.figsize']=(5,2)


        #An h5ad file is a hdf5 file with some additional 
        #structure  specifying  how to store AnnData 
        #objects.
        #What is an hdf5 file? 
        #It is a  Hierarchical Data Format 
        #(HDF) (HDF4, HDF5) designed to store and
        #organize large amounts of data.
        #Originally developed at the U.S.
        #National Center for Supercomputing Applications.

        self.results_file = 'write/pbmc3k.h5ad'

        #Generate annotated data file.
        self.source_path = 'data/filtered_gene_bc_matrices/hg19/'
        var_names = 'gene_symbols'
        self.adata = sc.read_10x_mtx(self.source_path,
                                var_names = var_names,
                                cache=True)
        print(self.adata)

    def check_for_repetitions(self):
        #Check if we have repetitions in the barcodes
        #or in the genes.

        #print(self.adata.obs)
        #print(self.adata.var)

        #This method can induce uniqueness.
        #self.adata.var_names_make_unique()

        print('>>>Checking for repetitions in:')

        print('Gene ids')
        if self.adata.var['gene_ids'].value_counts().gt(1).any():
            raise ValueError('Gene ids are not unique.')

        print('Barcodes')
        if self.adata.obs.index.value_counts().gt(1).any():
            raise ValueError('Barcodes are not unique.')

    def plot_genes_with_highest_exp(self, n=20):
        #First, we compute the total counts for each cell.
        #Then, we divide each entry by the total and 
        #convert that into a percentage.
        #Next, we plot the distribution of percentages for
        #those genes that had the largest medians.

        print('>>>Plotting the genes with the largest expression')
        save_path = '_top_' + str(n) + '_genes.png'
        #sc.pl.highest_expr_genes(self.adata,
                                 #n_top=n,
                                 #show= False,
                                 #save=save_path)

        #Let's do it by hand.
        #You can do it!
        total_counts_for_cells = self.adata.X.sum(axis=1)
        m = self.adata.X / total_counts_for_cells * 100

        df = pd.DataFrame(m.T, index = self.adata.var.index)
        medians = df.median(axis=1).sort_values(ascending=False)
        medians = medians.iloc[:n]
        #print(medians)
        df = df.loc[medians.index]

        factor = 1.5
        quantiles = df.quantile([0.25, 0.50, 0.75], axis=1)
        IQR = quantiles.iloc[-1] - quantiles.iloc[0]
        lower_bound = quantiles.iloc[0]  - factor * IQR
        upper_bound = quantiles.iloc[-1] + factor * IQR

        plt.close('all')
        fig,ax = plt.subplots()
        df = pd.melt(df.transpose(),
                     var_name = 'Gene',
                     value_name='% counts')

        sb.boxplot(data = df,
                   y='Gene',
                   x='% counts',
                   fliersize=1,
                   whis=factor,
                   ax = ax)

        #plt.ticklabel_format(style='sci',
                             #axis='both',
                             #scilimits=(0,0))
#
        #ax.set_xlim([0,2000])
        ax.set_title('Top expression in genes.')
        ax.set_ylabel('')
        save_path = os.path.join('.',
                                 'figures',
                                 'top_expression_in_genes.png')
        plt.tight_layout()
        fig.savefig(save_path)


    def apply_filter(self):
        pass
        print('>>>Filtering cells based on # of genes expressed')


    def keep_cells_with_at_least_n_genes(self,
                                         min_n_genes=200):
        txt = f'We require at least {min_n_genes} genes per cell.'
        print(txt)

        sc.pp.filter_cells(self.adata,
                           min_genes=min_n_genes,
                           inplace=True)


    def keep_genes_present_at_least_n_cells(self,
                                            min_n_cells=3):
        # We are going to filter the genes 
        # based on the minimum number of cells 
        # that expressed a given gene.

        txt  = '>>>Filtering genes based on # of cells'
        txt +=  ' having a given gene.'
        print(txt)

        min_n_cells = 3
        txt = f'We require at least {min_n_cells} cells per gene.'
        print(txt)
        sc.pp.filter_genes(self.adata,
                           min_cells=min_n_cells,
                           inplace=True)

    def find_mitochondrial_genes(self):

        mt_genes = self.adata.var_names.str.startswith('MT-')
        print('We have', sum(mt_genes), 'mitochondrial genes')

        #Add a column to the Genes data frame
        self.adata.var['mt'] = mt_genes
        mt_slice = self.adata.X[:, mt_genes]
        self.n_of_mt_reads = int(mt_slice.sum().sum())
        txt = 'We have {:d} reads associated to MT genes.'
        s = txt.format(self.n_of_mt_reads)
        print(s)

    def compute_qc_metrics(self):
        print('Computing QC metrics.')
        sc.pp.calculate_qc_metrics(self.adata,
                                   qc_vars=['mt'],
                                   percent_top=None,
                                   log1p=False,
                                   inplace=True)

    def plot_qc_metrics(self):
        print('Plotting QC metrics.')
        list_of_plots = ['n_genes_by_counts',
                         'total_counts',
                         'pct_counts_mt']

        save_path = '_QC_multi_panel.png'

        sc.pl.violin(self.adata, list_of_plots,
                     jitter=0.4, 
                     multi_panel=True, 
                     show=self.show_plot, 
                     save=save_path)


    def plot_count_dist_for_genes(self):

        print('Plot counts distribution for genes.')
        plt.close('all')
        fig,ax = plt.subplots()
        sb.histplot(self.adata.var,
                    x='total_counts',
                    bins=600,
                    kde=False,
                    ax = ax)

        plt.ticklabel_format(style='sci',
                             axis='both',
                             scilimits=(0,0))

        ax.set_xlim([0,2000])
        ax.set_title('Count dist. for genes')
        save_path = os.path.join('.',
                                 'figures',
                                 'count_dist_for_genes.png')
        fig.savefig(save_path)

    def plot_count_dist_for_cells(self):

        print('Plot counts distribution for cells.')
        plt.close('all')
        fig,ax = plt.subplots()
        sb.histplot(self.adata.obs,
                    x='total_counts',
                    bins=600,
                    kde=False,
                    ax = ax)

        #plt.ticklabel_format(style='sci',
                             #axis='both',
                             #scilimits=(0,0))

        #ax.set_xlim([0,2000])
        ax.set_title('Count dist. for cells')
        save_path = os.path.join('.',
                                 'figures',
                                 'count_dist_for_cells.png')
        fig.savefig(save_path)

    def plot_percent_mt_dist_for_cells(self):

        print('Plot percent of counts for mitochondrial genes distribution for cells.')
        plt.close('all')
        fig,ax = plt.subplots()
        sb.histplot(self.adata.obs,
                    x='pct_counts_mt',
                    bins=600,
                    kde=False,
                    ax = ax)

        #plt.ticklabel_format(style='sci',
                             #axis='both',
                             #scilimits=(0,0))

        #ax.set_xlim([0,2000])
        ax.set_title('% of counts (MT genes) for cells')
        save_path = os.path.join('.',
                                 'figures',
                                 'percent_mt_dist_for_cells.png')
        fig.savefig(save_path)

    def plot_var_vs_mean_for_genes(self):
        # We are going to plot, for each gene,
        # the variance against the mean.
        # We do that to illustrate that the associated
        # distribution is not Poisson.

        gene_column_means = self.adata.X.mean(axis=0)
        gene_column_means_sq = np.square(gene_column_means)
        sq_matrix = self.adata.X.power(2)
        gene_column_vars = sq_matrix.mean(axis=0) - gene_column_means_sq

        #print(gene_column_means.shape)
        #print(gene_column_vars.shape)

        x = np.asarray(gene_column_means).squeeze()
        y = np.asarray(gene_column_vars).squeeze()

        plt.close('all')
        fig,ax = plt.subplots()
        ax.scatter(x, y, s=5, c ='blue')
        ax.set_xlabel('Mean of counts for a given gene')
        ax.set_ylabel('Var of counts (gene)')
        plt.tight_layout()
        #ax.set_tite('Var vs mean of counts (gene)')
        save_path = './figures/variance_vs_mean_for_each_gene.png'
        fig.savefig(save_path)


    def run(self):

        self.check_for_repetitions()
        self.plot_genes_with_highest_exp()
        self.keep_cells_with_at_least_n_genes(200)
        self.keep_genes_present_at_least_n_cells(3)
        self.find_mitochondrial_genes()
        self.compute_qc_metrics()
        self.plot_count_dist_for_genes()
        self.plot_count_dist_for_cells()
        self.plot_percent_mt_dist_for_cells()
        self.plot_var_vs_mean_for_genes()

obj = Perseus()
obj.run()
