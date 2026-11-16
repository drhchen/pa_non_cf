source("results/baps_palette.R")
library(tidyverse); library(patchwork); library(ggrepel); library(ape); library(scales)
base <- normalizePath('results/3.6_LoF_analysis/results', mustWork=TRUE); out <- normalizePath('results/3.6_LoF_analysis/results/plotcase_figures', mustWork=FALSE)
dir.create(out,showWarnings=FALSE)
g <- read_tsv(file.path(base,'gene_LoF_summary.tsv'),show_col_types=FALSE)
e <- read_tsv(file.path(base,'LoF_events_with_metadata.tsv'),show_col_types=FALSE)
b <- read_tsv(file.path(base,'sample_LoF_burden.tsv'),show_col_types=FALSE)
focus <- c('mucA','mexB','lasR','oprD','mutS','mutL','ampD','algU','wspF','dacB')
# 1 frequency vs independent origins
p1d <- g %>% mutate(focus=gene_key %in% focus)
p1 <- ggplot(p1d,aes(samples,independent_origins,size=independent_ST,colour=independent_subclades)) + geom_point(alpha=.72) + geom_text_repel(data=filter(p1d,focus),aes(label=gene_key),size=3,box.padding=.4,point.padding=.3,segment.colour='grey45',max.overlaps=Inf) + scale_x_continuous(labels=comma) + scale_colour_gradientn(colours=unname(BAPS_COLORS), values=scales::rescale(1:8), breaks=1:8, labels=paste0("BAPS",1:8)) + scale_size_area(max_size=9) + labs(x='Isolates carrying a LoF variant',y='Independent LoF alleles',size='Independent STs',colour='Independent subclades',title='Clonal expansion versus recurrent LoF origins') + theme_classic(base_size=10)
# 2 mutation composition
m <- e %>% filter(gene_key %in% focus) %>% count(gene_key,mutation_class) %>% group_by(gene_key) %>% mutate(total=sum(n),prop=n/total) %>% ungroup() %>% mutate(gene_key=factor(gene_key,levels=rev(focus)))
p2 <- ggplot(m,aes(prop,gene_key,fill=mutation_class)) + geom_col(width=.72) + geom_text(data=filter(m,prop>.08),aes(label=n),position=position_stack(vjust=.5),size=3,colour='white') + scale_x_continuous(labels=percent,limits=c(0,1),expand=c(0,0)) + scale_fill_manual(values=c(frameshift='#2166AC',nonsense='#D6604D')) + labs(x='Proportion of LoF events',y=NULL,fill=NULL,title='Mutation-type composition of focal LoF genes') + theme_classic(base_size=10)
# 3 branch rank
br <- read_tsv(file.path(base,'focal_branch_normalized_permutation.tsv'),show_col_types=FALSE) %>% arrange(origin_rate_per_branch) %>% mutate(gene_key=factor(gene_key,levels=gene_key))
p3 <- ggplot(br,aes(origin_rate_per_branch,gene_key,size=samples,colour=independent_origins)) + geom_point() + scale_x_log10() + scale_colour_gradientn(colours=unname(BAPS_COLORS), values=scales::rescale(1:8), breaks=1:8, labels=paste0("BAPS",1:8)) + scale_size_area(max_size=9) + labs(x='Independent origins / terminal branch exposure (log scale)',y=NULL,size='Isolates',colour='Independent origins',title='Exploratory branch-exposure-normalized recurrence') + theme_classic(base_size=10)
fig <- (p1|p2)/p3 + plot_annotation(tag_levels='A',title='Supplementary analyses of recurrent LoF variation',theme=theme(plot.title=element_text(face='bold',size=15)))
ggsave(file.path(out,'Supp_Figure_LoF_recurrence_mechanisms_plotcase.pdf'),fig,width=12,height=10,device=cairo_pdf)
ggsave(file.path(out,'Supp_Figure_LoF_recurrence_mechanisms_plotcase.png'),fig,width=12,height=10,dpi=400)
# 4 tree-order heatmap, using tree tip order
tr <- read.tree('/Users/hc/Downloads/投稿文章/pa_non_cf/alldata/tree/panaroo-core-genome-tree/panaroo_483_core.snp.filter.aln.treefile')
ord <- tr$tip.label
mat <- expand_grid(sample=ord,gene_key=focus) %>% left_join(e %>% filter(gene_key %in% focus) %>% distinct(sample,gene_key,mutation_class),by=c('sample','gene_key')) %>% mutate(value=case_when(mutation_class=='frameshift'~'Frameshift',mutation_class=='nonsense'~'Nonsense',is.na(mutation_class)~'None',TRUE~'Both'))
mat$sample=factor(mat$sample,levels=rev(ord)); mat$gene_key=factor(mat$gene_key,levels=focus)
p4 <- ggplot(mat,aes(gene_key,sample,fill=value))+geom_tile() + scale_fill_manual(values=c(None='white',Frameshift='#2166AC',Nonsense='#D6604D',Both='#984EA3'),drop=FALSE) + labs(x='Focal LoF gene',y='Core-genome tree tip order',fill=NULL,title='Focal LoF variants in phylogenetic tip order') + theme_minimal(base_size=9)+theme(axis.text.y=element_blank(),panel.grid=element_blank(),legend.position='right')
ggsave(file.path(out,'Supp_Figure_LoF_phylogeny_ordered_heatmap_plotcase.pdf'),p4,width=7,height=11,device=cairo_pdf)
ggsave(file.path(out,'Supp_Figure_LoF_phylogeny_ordered_heatmap_plotcase.png'),p4,width=7,height=11,dpi=400)
