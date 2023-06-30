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

        #Figure path
        sc.settings.figdir = os.path.join('.', 'figures')


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

        self.results_file = '/home/javier/Documents/perseus/write/pbmc3k.h5ad'

        #Generate annotated data file.
        self.source_path = '/home/javier/Documents/perseus/data/filtered_gene_bc_matrices/hg19'
        var_names = 'gene_symbols'
        self.adata = sc.read_10x_mtx(self.source_path,
                                var_names = var_names,
                                cache=True)


        print(self.adata)

    def compute_total_counts(self):
        '''
        Calculate the total counts for each cell.
        Calculate the global total count, which is the
        sum of the entries of the total counts for each cell.
        '''

        self.total_counts_for_cells = self.adata.X.sum(axis=1)
        self.global_total_counts = self.total_counts_for_cells.sum()

    def check_for_repetitions(self):
        '''Check if we have repetitions in the barcodes
        or in the genes.
        '''

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
        ''' First, we compute the total counts for each cell.
        Then, we divide each entry by the total and 
        convert that into a percentage.
        Next, we plot the distribution of percentages for
        those genes that had the largest medians.

        param int n: Number of genes to show.

        '''

        self.compute_total_counts()

        print('>>>Plotting the genes with the largest expression')
        save_path = '_top_' + str(n) + '_genes.png'
        #sc.pl.highest_expr_genes(self.adata,
                                 #n_top=n,
                                 #show= False,
                                 #save=save_path)

        #Let's do it by hand.
        #You can do it!
        m = self.adata.X / self.total_counts_for_cells * 100

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

        ax.set_title('Top expression in genes.')
        ax.set_ylabel('')
        save_path = os.path.join('.',
                                 'figures',
                                 'top_expression_in_genes.png')
        plt.tight_layout()
        fig.savefig(save_path)


    def keep_cells_with_at_least_n_genes(self,
                                         min_n_genes=200):
        '''Eliminate cells that have less than n genes.

        :param int min_n_genes: Minimum number of genes in a cell.
        
        '''
        txt = f'We require at least {min_n_genes} genes per cell.'
        print(txt)

        sc.pp.filter_cells(self.adata,
                           min_genes=min_n_genes,
                           inplace=True)


    def keep_genes_present_in_at_least_n_cells(self,
                                            min_n_cells=3):
        '''Filter the genes based on the minimum number of cells 
        that expressed a given gene.

        :param int min_n_cells: Minimum number of cells
        '''

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
        '''
        Create a column in the Genes data frame specifying
        which genes are mitochondrial.
        Also print how many of the total counts correspond to
        mitochondrial genes.
        '''

        mt_genes = self.adata.var_names.str.startswith('MT-')
        print('We have', sum(mt_genes), 'mitochondrial genes')

        #Add a column to the Genes data frame
        self.adata.var['mt'] = mt_genes
        mt_slice = self.adata.X[:, mt_genes]
        self.global_mt_reads = int(mt_slice.sum().sum())
        self.global_pct_mt_reads = self.global_mt_reads / self.global_total_counts * 100
        txt = 'We have {:d} reads associated to MT genes.'
        s = txt.format(self.global_mt_reads)
        print(s)

        txt = 'We have {:2f} % of reads associated to MT genes.'
        s = txt.format(self.global_pct_mt_reads)
        print(s)

    def compute_qc_metrics(self):
        '''
        Compute quality control metrics for cells and genes.
        These metrics are stored in the dataframes obs=cells
        and var=genes.
        Mitochondrial information is also included but first 
        we need to call :py:func:`find_mitochondrial_genes`.
        '''
        print('Computing QC metrics.')
        sc.pp.calculate_qc_metrics(self.adata,
                                   qc_vars=['mt'],
                                   percent_top=None,
                                   log1p=False,
                                   inplace=True)

    def plot_standard_qc_metrics(self):
        '''
        Generate 3 plots showing the following:
        #. For each cell (dot), the number of genes that have a
        positive count.
        #. The counts for each cell.
        #. The % of counts that correspond to mitochondrial
        genes for each cell.
        '''
        print('Plotting QC metrics.')
        list_of_plots = ['n_genes_by_counts',
                         'total_counts',
                         'pct_counts_mt']

        save_path = '_QC_multi_panel.png'

        #print(mpl.rcParams['figure.figsize'])
        #mpl.rcParams['figure.figsize']=(5,2)
        sc.pl.violin(self.adata, list_of_plots,
                     jitter=0.4, 
                     multi_panel=True, 
                     show=self.show_plot, 
                     save=save_path)


    def plot_count_dist_for_genes(self):
        '''
        Show the count distribution for each gene.
        '''

        print('Plot counts distribution for genes.')
        plt.close('all')
        fig,ax = plt.subplots()
        sb.histplot(self.adata.var,
                    x='total_counts',
                    #bins=600,
                    kde=False,
                    ax = ax)

        plt.ticklabel_format(style='sci',
                             axis='both',
                             scilimits=(0,0))

        #alert
        ax.set_xlim([0,2000])
        ax.set_ylabel('Frequency')
        ax.set_title('Count dist. for genes')
        save_path = os.path.join('.',
                                 'figures',
                                 'count_dist_for_genes.png')
        fig.savefig(save_path)

    def global_n_cells(self):
        '''
        Provide the current number of active cells.
        '''
        return self.adata.obs.shape[0]

    def plot_gene_dist_for_cells(self):
        '''
        We show the distribution of the number of genes
        with positive counts for each cell.
        '''

        print('Plotting gene # distribution for cells.')
        plt.close('all')
        fig,ax = plt.subplots()
        sb.histplot(self.adata.obs,
                    x='n_genes_by_counts',
                    bins=600,
                    kde=False,
                    ax = ax)

        #plt.ticklabel_format(style='sci',
                             #axis='both',
                             #scilimits=(0,0))

        #ax.set_xlim([0,2000])
        ax.set_ylabel('Frequency')
        ax.set_title('Gene # dist. for cells')
        save_path = os.path.join('.',
                                 'figures',
                                 'gene_n_dist_for_cells.png')
        fig.savefig(save_path)

    def plot_count_dist_for_cells(self):
        '''
        Show the count distribution for each cell.
        '''

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
        ax.set_ylabel('Frequency')
        ax.set_title('Count dist. for cells')
        save_path = os.path.join('.',
                                 'figures',
                                 'count_dist_for_cells.png')
        fig.savefig(save_path)

    def plot_percent_mt_dist_for_cells(self):
        '''Plot the percent of counts corresponding to 
        mitochondrial genes for each cell as a distribution.
        '''

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
        ax.set_ylabel('Frequency')
        ax.set_title('% of counts (MT genes) for cells')
        save_path = os.path.join('.',
                                 'figures',
                                 'percent_mt_dist_for_cells.png')
        fig.savefig(save_path)

    def plot_var_vs_mean_for_genes(self):
        '''Plot for each gene,
         the variance against the mean for the counts.
         We do that to illustrate that the associated
         distribution is not Poisson. We say that
         :math:`X\\sim\\textrm{Poisson}(\\lambda)` if
         :math:`P(X=k) = \\dfrac{e^{-\\lambda} \\lambda^k}{k!}`
        '''

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

    def plot_genes_vs_total_counts_vs_pct_mt(self):
        '''
        The y-axis indicates the number of genes with a
        positive count for a given cell. The x-axis is
        the number of counts for a given cell.
        The coloring is a function of the percentage of
        mitochondrial counts for a given cell.
        '''

        print('Plotting genes vs counts vs pct mt')
        save_path = '_genes_counts_mt_color.png'
        sc.pl.scatter(self.adata,
                      x='total_counts',
                      y='n_genes_by_counts',
                      title='# of genes vs # count vs % mt',
                      color='pct_counts_mt',
                      save=save_path,
                      show=self.show_plot)

    def count_IQR_filtration_for_cells(self, factor=1.5):
        '''
        Filtrate cells that have counts outside the range
        :math:`(Q(25\\%)-f\\cdot IQR, Q(75\\%)+f\\cdot IQR)`
        where :math:`f` is a range-elongating factor.
        '''
        print('Filtrating cells by counts')
        q = self.adata.obs.total_counts.quantile([0.25, 0.75])
        IQR = q.iloc[-1] - q.iloc[0]
        upper = q.iloc[-1] + factor * IQR
        lower = q.iloc[0]  - factor * IQR
        mask = upper < self.adata.obs.total_counts
        mask |= self.adata.obs.total_counts < lower
        f = mask.sum()
        f_pct = f / self.global_n_cells() * 100
        print(f'There are {f} cells outside the range.')
        print(f'This corresponds to {f_pct:.2f} % of all cells.')
        print('Removing the outliers.')
        self.adata = self.adata[~mask, :]


    def pct_mt_IQR_filtration_for_cells(self, factor=1.5):
        '''
        Filtrate cells that have a % of counts corresponding
        to mitochondrial genes outside the range
        :math:`(Q(25\\%)-f\\cdot IQR, Q(75\\%)+f\\cdot IQR)`
        where :math:`f` is a range-elongating factor.
        '''
        print('Filtrating cells by % counts (mitochondrial genes)')
        q = self.adata.obs.pct_counts_mt.quantile([0.25, 0.75])
        IQR = q.iloc[-1] - q.iloc[0]
        upper = q.iloc[-1] + factor * IQR
        lower = q.iloc[0]  - factor * IQR
        mask = upper < self.adata.obs.pct_counts_mt
        mask |= self.adata.obs.pct_counts_mt < lower
        f = mask.sum()
        f_pct = f / self.global_n_cells() * 100
        print(f'There are {f} cells outside the range.')
        print(f'This corresponds to {f_pct:.2f} % of all cells.')
        self.adata = self.adata[~mask, :]

    def gene_IQR_filtration_for_cells(self, factor=1.5):
        '''
        Filtrate cells that expressed a number of genes
        outside the range
        :math:`(Q(25\\%)-f\\cdot IQR, Q(75\\%)+f\\cdot IQR)`
        where :math:`f` is a range-elongating factor.
        '''
        print('Filtrating cells by # of genes')
        q = self.adata.obs.n_genes_by_counts.quantile([0.25, 0.75])
        IQR = q.iloc[-1] - q.iloc[0]
        upper = q.iloc[-1] + factor * IQR
        lower = q.iloc[0]  - factor * IQR
        mask = upper < self.adata.obs.n_genes_by_counts
        mask |= self.adata.obs.n_genes_by_counts < lower
        f = mask.sum()
        f_pct = f / self.global_n_cells() * 100
        print(f'There are {f} cells outside the range.')
        print(f'This corresponds to {f_pct:.2f} % of all cells.')
        self.adata = self.adata[~mask, :]

    def normalize_by_medians(self):
        '''
        Compute the median count for each row(=cell)
        and divide the corresponding row by the
        median.
        '''
        print(self.adata.X)
        medians = self.adata.X.median(axis=1)
        self.adata.X = self.adata.X / medians
        print(self.adata.X)

    def compute_z_scores_using_cells(self):
        '''
        First we compute :math:`\\mu` and
        :math:`\\sigma`, which are the
        mean and variance of counts for each row(=cell)
        Then, for each entry, we compute the z-score.
        Recall that :math:`\\sigma^2 = E(X^2) - E^2(X)`.
        We use this fact to compute the variance.
        '''
        sq_matrix = self.adata.X.power(2)
        mu_X    = self.adata.X.mean(axis=1)
        sq_mu_X = np.multiply(mu_X,mu_X)
        mu_X_sq = sq_matrix.mean(axis=1)
        sigma_sq = mu_X_sq - sq_mu_X
        sigma = np.sqrt(sigma_sq)
        indices = self.adata.X.nonzero()
        mu_X = np.asarray(mu_X).squeeze()
        self.adata.X[indices] -= mu_X[indices[0]]
        self.adata.X /= sigma

    def run(self):
        '''
        Execute the pipeline to process the count matrix.
        #. Filtration
        #. Normalization
        #. Dimensionality reduction
        #. Visualization
        '''

        self.check_for_repetitions()
        self.plot_genes_with_highest_exp()
        self.keep_cells_with_at_least_n_genes(200)
        self.keep_genes_present_in_at_least_n_cells(3)
        self.find_mitochondrial_genes()
        self.compute_qc_metrics()
        self.plot_standard_qc_metrics()
        #self.plot_count_dist_for_genes()
        #self.plot_count_dist_for_cells()
        #self.plot_gene_dist_for_cells()
        #self.plot_percent_mt_dist_for_cells()
        #self.plot_var_vs_mean_for_genes()
        self.plot_genes_vs_total_counts_vs_pct_mt()
        self.gene_IQR_filtration_for_cells()
        self.count_IQR_filtration_for_cells()
        self.pct_mt_IQR_filtration_for_cells()
        #self.normalize_by_medians()
        self.compute_z_scores_using_cells()

