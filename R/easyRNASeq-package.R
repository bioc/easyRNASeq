# to re-build the documentation and namespace
# library(roxygen2)
# roxygenize("../easyrnaseq-devel/",roclets=c('rd', 'collate', 'namespace'),clean=TRUE)

# to update the package versions
# pkg <- c("Biobase","BiocGenerics","BiocParallel","biomaRt","Biostrings",
#          "DESeq","edgeR","GenomeInfoDb","genomeIntervals",
#          "GenomicAlignments","GenomicRanges","SummarizedExperiment",
#          "IRanges","LSD","Rsamtools","S4Vectors","ShortRead",
#          "BiocStyle","rappdirs,
#          "BSgenome",
#          "BSgenome.Dmelanogaster.UCSC.dm3","GenomicFeatures",
#          "RUnit")
# pkg[! pkg %in% rownames(installed.packages())]
# installed.packages()[pkg,"Version"]

##==========================
## package details
##==========================
#' Count summarization and normalization pipeline for Next Generation Sequencing data.
#'
#' Offers functionalities to summarize read counts per feature of interest, e.g. exons, transcripts, genes, etc.
#' Offers functionalities to normalize the summarized counts using 3rd party packages like \code{\link[DESeq:newCountDataSet]{DESeq}}
#' or \code{\link[edgeR:DGEList]{edgeR}}.
#'
#' @section Methods:
#' The main function \code{\link[easyRNASeq:easyRNASeq-easyRNASeq]{easyRNASeq}} will summarize the counts per
#' feature of interest, for as many samples as provided and will return a
#' count matrix (N*M) where N are the features and M the samples.
#' This data can be corrected to \pkg{RPKM} in which case
#' a matrix of corrected value is returned instead, with the same dimensions.
#' Alternatively a \code{\linkS4class{RangedSummarizedExperiment}} can be returned and this
#' is expected to be the default in the upcoming version of easyRNASeq (as of 1.5.x).
#' If the necessary sample
#' information are provided, the data can be normalized using either \code{\link[DESeq:newCountDataSet]{DESeq}}
#' or \code{\link[edgeR:DGEList]{edgeR}} and the corresponding package object returned.
#' For more insider details, and step by step functions, see:
#' \tabular{ll}{
#' 	\code{\link[easyRNASeq:ShortRead-methods]{ShortRead methods}} for pre-processing the data.
#' 	\code{\link[easyRNASeq:easyRNASeq-annotation-methods]{easyRNASeq annotation methods}} for getting the annotation.
#' 	\code{\link[easyRNASeq:easyRNASeq-coverage-methods]{easyRNASeq coverage methods}} for computing the coverage from a Short Read Alignment file.
#' 	\code{\link[easyRNASeq:easyRNASeq-summarization-methods]{easyRNASeq summarization methods}} for summarizing the data.
#' 	\code{\link[easyRNASeq:easyRNASeq-correction-methods]{easyRNASeq correction methods}} for correcting the data (i.e. generating RPKM).
#' 	\code{\link[easyRNASeq:edgeR-methods]{edgeR methods}} for post-processing the data.
#' 	\code{\link[easyRNASeq:DESeq-methods]{DESeq methods}} for post-processing the data.
#' 	}
#'
#' @name easyRNASeq package
#' @rdname easyRNASeq-package
#' @aliases easyRNASeq-package assay type BamFileList BamFileList-class IRanges
#' GRanges-class GRanges SRFilterResult SummarizedExperiment RangedSummarizedExperiment-class
#' chromosomeFilter compose nFilter
#' @docType package
#' @author Nicolas Delhomme, Bastian Schiffthaler, Ismael Padioleau
#' @keywords package
#' @seealso
#' 	The class RNAseq specification:
#' 	\code{\linkS4class{RNAseq}}
#'
#'   The default output class specification:
#'   \code{\linkS4class{RangedSummarizedExperiment}}
#'
#' 	The imported packages:
#' 	\code{\link[biomaRt:useMart]{biomaRt}}
#'   \code{\link[BiocParallel:BiocParallel-package]{BiocParallel}}
#' 	\code{\link[edgeR:DGEList]{edgeR}}
#' 	\code{\link[genomeIntervals:Genome_intervals_stranded-class]{genomeIntervals}}
#' 	\code{\link[Biostrings:XString-class]{Biostrings}}
#' 	\code{\link[BSgenome:BSgenome-class]{BSgenome}}
#' 	\code{\link[DESeq:newCountDataSet]{DESeq}}
#' 	\code{\link[GenomicRanges:GRanges-class]{GenomicRanges}}
#' 	\code{\link[IRanges:IRanges-constructor]{IRanges}}
#' 	\code{\link[Rsamtools:scanBam]{Rsamtools}}
#' 	\code{\link[ShortRead:readAligned]{ShortRead}}
#'
#' 	The suggested packages:
#' 	\code{\link[parallel:makeCluster]{parallel}}
#' 	\code{\link[GenomicFeatures:TxDb-class]{GenomicFeatures}}
#'
#'   The following classes and functions that are made available from
#'   other packages:
#'   \itemize{
#'     \item{Classes}{
#'       \code{\link[Rsamtools:BamFile-class]{BamFileList-class}}
#'       \code{\linkS4class{CountDataSet}}
#'       \code{\linkS4class{RangedSummarizedExperiment}}
#'     }
#'     \item{Functions/Methods}{
#'       \code{\link[DESeq:estimateDispersions]{DESeq estimate size factor and
#'       estimate dispersion functions}}
#'       \code{\link[SummarizedExperiment:RangedSummarizedExperiment-class]{
#'         The RangedSummarizedExperiment assay accessor}
#'       }
#'       The locfit function \code{\link[locfit]{locfit}}
#'       The BamFileList constructor \code{\link[Rsamtools:BamFile-class]{BamFileList-class}}
#'       The IRanges constructor \code{\link[IRanges]{IRanges-constructor}}
#'       For the SRFilterResult,
#'       chromosomeFilter, compose and nFilter methods\code{\link[ShortRead]{srFilter}}
#'     }
#'   }
#'
#' @examples
#'  # the data
#'  tdir <- tutorialData()
#'
#'   # get the example annotation file - we retrieve a gtf file from GitHub
#'   annot <- fetchData("Drosophila_melanogaster.BDGP5.77.with-chr.gtf.gz")
#'
#'   # create the AnnotParam
#'   annotParam <- AnnotParam(
#'     datasource=annot,
#'     type="gtf")
#'
#'   # create the synthetic transcripts
#'   annotParam <- createSyntheticTranscripts(annotParam,verbose=FALSE)
#'
#'  # create the RnaSeqParam
#'  rnaSeqParam <- RnaSeqParam(annotParam=annotParam,countBy="gene")
#'
#'   # get the bamfiles (from the Bioc cache in this example)
#'  filenames <- dir(tdir,pattern="[A,T].*\\.bam$",full.names=TRUE)
#'  indexnames <- sapply(paste0(sub(".*_","",basename(filenames)),".bai"),fetchData)
#'  bamFiles <- getBamFileList(filenames,indexnames)
#'
#'  # get a RangedSummarizedExperiment containing the counts table
#'   sexp <- simpleRNASeq(
#'       bamFiles=bamFiles,
#'       param=rnaSeqParam,
#'       verbose=TRUE
#'   )
#'
#'   # get the counts
#'   assays(sexp)$genes
#'
#'
NULL

##==========================
# To define the NAMESPACE
##==========================
# library(codetoolsBioC)
# library(easyRNASeq)
# writeNamespaceImports("easyRNASeq")
# import classes
#' @importClassesFrom Biostrings DNAStringSet
#' @importClassesFrom DESeq CountDataSet
#' @importClassesFrom edgeR DGEList
#' @importClassesFrom genomeIntervals Genome_intervals Genome_intervals_stranded
#' @importClassesFrom GenomicAlignments GAlignments GAlignmentPairs
#' @importClassesFrom GenomicRanges GRanges GRangesList
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
#' @importClassesFrom IRanges RleList
#' @importClassesFrom S4Vectors Annotated Vector DataFrame SimpleList
#' @importClassesFrom methods ANY character "function" integer
#' list matrix missing numeric vector
#' @importClassesFrom Rsamtools BamFile BamFileList
#' @importClassesFrom ShortRead AlignedRead ShortReadQ
# import S4 methods
#' @importMethodsFrom Biobase fData varMetadata
#' @importFrom BiocFileCache bfcadd bfcdownload bfcneedsupdate bfcquery bfcrpath
#' @importMethodsFrom BiocGenerics annotation cbind clusterApply
#' colnames counts duplicated estimateDispersions estimateSizeFactors
#' eval fileName get intersect lapply match order path paste pmax rbind
#' rownames sapply strand "strand<-" table unique
#' @importMethodsFrom Biostrings type
#' @importMethodsFrom DESeq estimateSizeFactors estimateDispersions
#' @importMethodsFrom genomeIntervals readGff3 writeGff3
#' @importMethodsFrom GenomeInfoDb seqinfo seqlengths "seqlengths<-"
#' seqlevels "seqlevels<-" seqnames "seqnames<-"
#' @importMethodsFrom GenomicAlignments cigar summarizeOverlaps
#' @importMethodsFrom GenomicRanges grglist
#' @importMethodsFrom IRanges as.list as.matrix
#' "colnames<-" countOverlaps coverage end "end<-" findOverlaps
#' gsub mean median narrow nchar ranges reduce "rownames<-" space
#' start "start<-" sub  tolower unlist values which width
#' @importMethodsFrom methods coerce initialize show
#' @importMethodsFrom Rsamtools asMates "asMates<-" countBam scanBam
#' scanBamHeader ScanBamParam yieldSize "yieldSize<-"
#' @importMethodsFrom S4Vectors "%in%" as.table elementMetadata
#' "elementMetadata<-" elementNROWS levels mcols metadata
#' "metadata<-" Rle runLength runsum runValue split substr
#' @importMethodsFrom ShortRead chromosome id position readAligned
#' srdistance sread srFilter writeFastq
#' @importMethodsFrom SummarizedExperiment assay assays "assays<-"
#' colData "colData<-" rowRanges "rowRanges<-"
# import methods
#' @importFrom biomaRt getBM listDatasets useDataset useMart
#' @importFrom BiocParallel MulticoreParam SerialParam
#' @importFrom Biostrings DNAStringSet
#' @importFrom DESeq fitInfo newCountDataSet
#' @importFrom edgeR calcNormFactors DGEList estimateCommonDisp
#' estimateTagwiseDisp maPlot plotMDS.DGEList plotMeanVar plotBCV
#' @importFrom genomeIntervals getGffAttribute parseGffAttributes
#' @importFrom GenomicAlignments GAlignments readGAlignments readGAlignmentPairs
#' @importFrom GenomicRanges GRanges GRangesList
#' @importFrom graphics abline axis axTicks boxplot grid hist legend lines
#' mtext par plot rect
#' @importFrom IRanges IRanges IRangesList LogicalList
#' SplitDataFrameList RleList
#' @importFrom locfit locfit lp
#' @importFrom LSD heatscatter
#' @importFrom methods as extends is new
#' @importFrom parallel makePSOCKcluster parLapply stopCluster
#' @importFrom rappdirs user_cache_dir
#' @importFrom Rsamtools BamFileList bamFlagTest index scanBamFlag
#' @importFrom S4Vectors endoapply DataFrame queryHits SimpleList
#' @importFrom ShortRead alignData chromosomeFilter compose nFilter
#' SRFilterResult
#' @importFrom stats aggregate na.omit
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom utils combn str packageVersion
# and export!
#' @exportClass BamFileList GRanges RangedSummarizedExperiment
#' @exportMethod assay assays colData estimateDispersions estimateSizeFactors fileName metadata rowRanges seqinfo seqlengths seqlevels "seqlevels<-" seqnames "seqnames<-" split srFilter width writeFastq writeGff3
#' @export alignData basename chromosomeFilter compose BamFileList GRanges IRanges locfit lp newCountDataSet nFilter readAligned SRFilterResult SummarizedExperiment
NULL

##==========================
# To detail the deprecation
##==========================
# nothing at the time - copy the defunct block if needed and
# replace defunct by deprecated
##==========================
# To detail defunct function
##==========================
#' The following function are defunct:
#' \itemize{
#' \item \code{\link[easyRNASeq:easyRNASeq-easyRNASeq]{easyRNASeq}}
#' \item \code{\link[easyRNASeq:easyRNASeq-coverage-methods]{fetchCoverage}}
#' \item \code{fetchAnnotation}
#' \item \code{knownOrganisms}
#' \item \code{plotDispersionEstimates,DGEList-method}
#' }
#'
#' \itemize{
#' \item The \code{plotDispersionEstimates,DGEList-method}
#' function is superseded by the \code{\link[edgeR:plotBCV]{plotBCV}} function
#' as the \pkg{edgeR} DGEList object structure changed
#' \item The \code{\link[easyRNASeq:easyRNASeq-easyRNASeq]{easyRNASeq}} function is superseded by the
#' \code{\link[easyRNASeq:easyRNASeq-simpleRNASeq]{simpleRNASeq}} function to consolidate and
#' prune the overall package. The changes are based on user comments and on the
#' general standardization occuring in the field.
#' \item The \code{\link[easyRNASeq:easyRNASeq-coverage-methods]{fetchCoverage}} function only had two
#' parameters deprecated as the consequence of the package consolidation. As the scanBam
#' function is not called directly anymore but through higher level functions (from the
#' GenomicRanges package), the 'what' and 'isUnmappedQuery' parameters were obsolete.
#' }
#' @aliases plotDispersionEstimates,DGEList-method fetchAnnotation knownOrganisms
#' easyRNASeq easyRNASeq,RNAseq-method fetchCoverage,RNAseq-method fetchCoverage
#' organismName<- organismName organismName,RNAseq-method organismName<-,RNAseq-method
#' @name Defunct functions
#' @rdname easyRNASeq-defunct
NULL

##==========================
# To detail dataset
##==========================
#' Dataset included in the package
#'
#' The package contains a dataset from the \emph{Robinson, Delhomme et al., 2014}
#' publication.
#' \itemize{
#' \item{RobinsonDelhomme2014}{a normalised expression count table. This dataset was
#' generated from 17 \emph{Populus tremula} - Eurasian
#' aspen - trees used to assess the sexual dimorphism of this dioecious species. This
#' count matrix has been generating following published pre-processing
#' guidelines - see \url{http://www.epigenesys.eu/en/protocols/bio-informatics/1283-guidelines-for-rna-seq-data-analysis} -
#' and the resulting HTSeq files have been collated and the obtained raw count
#' matrix submitted to a variance stabilising transformation. Subsequently, the
#' values have been transformed so that the minimal vst values - that
#' corresponds to an absence of expression - is 0. Hence the counts in the
#' matrix are library-size normalized, variance stabilised expression values,
#' with a minimal value of 0.}
#' }
#' @aliases easyRNASeq-datasets
#' RobinsonDelhomme2014
#' @name easyRNASeq-datasets
#' @rdname easyRNASeq-datasets
#' @keywords data
NULL

##==========================
# To detail hardcoded parameters
##==========================
#' Objects created when the package is attached.
#'
#' The package creates the following objects when attached
#' \itemize{
#' \item GTF.FIELDS
#' \item ANNOTATION.TYPE
#' \item TUTORIAL.DATA
#' \item VIGNETTE.DATA
#' }
#'
#' These objects hold the following information
#' \itemize{
#' \item GTF.FIELDS \code{c("gene_id","transcript_id","exon_id","gene_name")}
#' \item ANNOTATION.TYPE \code{c(mRNA="mRNA",exon="exon")}
#' \item TUTORIAL.DATA: The list of files needed for the help and test
#' pages
#' \item VIGNETTE.DATA: The list of files needed for the vignette
#' }
#' and are designed as global variables to expose the
#' fact that they are hardcoded. These exist as
#' placeholder in case a user would require different
#' values for these.
#'
#' @aliases easyRNASeq-global-variables .onAttach ANNOTATION.TYPE GTF.FIELDS
#' TUTORIAL.DATA VIGNETTE.DATA
#' @name easyRNASeq-global-variables
#' @rdname easyRNASeq-global-variables
#' @param libname a character string giving the library directory where the package defining the namespace was found.
#' @param pkgname a character string giving the name of the package.
#' @seealso \code{\link[base:ns-hooks]{.onAttach}} in the \code{base} package.
#' @keywords internal
globalVariables("GTF.FIELDS")
globalVariables("ANNOTATION.TYPE")
globalVariables("TUTORIAL.DATA")
globalVariables("VIGNETTE.DATA")
".onAttach" <- function(libname,pkgname){
  assign("GTF.FIELDS",c("gene_id","transcript_id","exon_id",
                        "gene_name"),
         envir=as.environment("package:easyRNASeq"))
  assign("ANNOTATION.TYPE",c(mRNA="mRNA",exon="exon"),
         envir=as.environment("package:easyRNASeq"))
  .gitHubURL <- "https://github.com/UPSCb/UPSCb/raw/master/tutorial/easyRNASeq"
  assign("TUTORIAL.DATA",
         c(ACACTG.bam=file.path(.gitHubURL,"ACACTG.bam"),
           ACTAGC.bam=file.path(.gitHubURL,"ACTAGC.bam"),
           ATGGCT.bam=file.path(.gitHubURL,"ATGGCT.bam"),
           TTGCGA.bam=file.path(.gitHubURL,"TTGCGA.bam"),
           ACACTG.bam.bai=file.path(.gitHubURL,"ACACTG.bam.bai"),
           ACTAGC.bam.bai=file.path(.gitHubURL,"ACTAGC.bam.bai"),
           ATGGCT.bam.bai=file.path(.gitHubURL,"ATGGCT.bam.bai"),
           TTGCGA.bam.bai=file.path(.gitHubURL,"TTGCGA.bam.bai"),
           "Drosophila_melanogaster.BDGP5.77.with-chr.gtf.gz"=
               file.path(.gitHubURL,"Drosophila_melanogaster.BDGP5.77.with-chr.gtf.gz"),
           gAnnot.rda=file.path(.gitHubURL,"gAnnot.rda"),
           "Dmel-mRNA-exon-r5.52.gff3.gz"=file.path(.gitHubURL,"Dmel-mRNA-exon-r5.52.gff3.gz")),
  envir=as.environment("package:easyRNASeq"))
  assign("VIGNETTE.DATA",
         c("Ptrichocarpa_210_v3.0_gene_exons.gff3.gz"=
             file.path(.gitHubURL,"Ptrichocarpa_210_v3.0_gene_exons.gff3.gz"),
           "md5.txt"=file.path(.gitHubURL,"md5.txt"),
           "202_subset_sortmerna_trimmomatic_sorted.bam"=
             file.path(.gitHubURL,"202_subset_sortmerna_trimmomatic_sorted.bam"),
           "207_subset_sortmerna_trimmomatic_sorted.bam"=
             file.path(.gitHubURL,"207_subset_sortmerna_trimmomatic_sorted.bam"),
           "213.1_subset_sortmerna_trimmomatic_sorted.bam"=
             file.path(.gitHubURL,"213.1_subset_sortmerna_trimmomatic_sorted.bam"),
           "221_subset_sortmerna_trimmomatic_sorted.bam"=
             file.path(.gitHubURL,"221_subset_sortmerna_trimmomatic_sorted.bam"),
           "226.1_subset_sortmerna_trimmomatic_sorted.bam"=
             file.path(.gitHubURL,"226.1_subset_sortmerna_trimmomatic_sorted.bam"),
           "229.1_subset_sortmerna_trimmomatic_sorted.bam"=
             file.path(.gitHubURL,"229.1_subset_sortmerna_trimmomatic_sorted.bam"),
           "202_subset_sortmerna_trimmomatic_sorted.bam.bai"=
             file.path(.gitHubURL,"202_subset_sortmerna_trimmomatic_sorted.bam.bai"),
           "207_subset_sortmerna_trimmomatic_sorted.bam.bai"=
             file.path(.gitHubURL,"207_subset_sortmerna_trimmomatic_sorted.bam.bai"),
           "213.1_subset_sortmerna_trimmomatic_sorted.bam.bai"=
             file.path(.gitHubURL,"213.1_subset_sortmerna_trimmomatic_sorted.bam.bai"),
           "221_subset_sortmerna_trimmomatic_sorted.bam.bai"=
             file.path(.gitHubURL,"221_subset_sortmerna_trimmomatic_sorted.bam.bai"),
           "226.1_subset_sortmerna_trimmomatic_sorted.bam.bai"=
             file.path(.gitHubURL,"226.1_subset_sortmerna_trimmomatic_sorted.bam.bai"),
           "229.1_subset_sortmerna_trimmomatic_sorted.bam.bai"=
             file.path(.gitHubURL,"229.1_subset_sortmerna_trimmomatic_sorted.bam.bai")),
         envir=as.environment("package:easyRNASeq"))
}
