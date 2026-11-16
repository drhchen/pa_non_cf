source("results/baps_palette.R")
BAPS_PLOT_COLORS <- BAPS_COLORS[paste0("BAPS", c(1,2,3,4,5,6,7,8))]
suppressPackageStartupMessages({library(tidyverse);library(patchwork);library(igraph);library(ggraph)})
root <- 'results/3.5_mobilome/results_fastbaps'; out <- file.path(root,'figures_final'); dir.create(out,recursive=TRUE,showWarnings=FALSE)
ord <- c('BAPS1','BAPS6','BAPS7','BAPS8')
inv <- read.delim(file.path(root,'element_inventory_final.tsv'),check.names=FALSE) %>% mutate(BAPS=factor(BAPS,levels=ord))
# A: only observed elements; sum all elements within sample x type; no zero imputation
len <- inv %>% filter(BAPS %in% ord, element_type != 'compositional_outlier') %>% group_by(sample_id,BAPS,element_type) %>% summarise(total_length_kb=sum(length_bp,na.rm=TRUE)/1000,.groups='drop')
pA <- ggplot(len,aes(BAPS,total_length_kb,colour=BAPS)) + geom_boxplot(outlier.shape=NA,width=.62,alpha=.12) + geom_jitter(width=.12,height=0,size=.7,alpha=.55) + facet_wrap(~element_type,scales='free_y',ncol=5) + scale_colour_manual(values=BAPS_PLOT_COLORS, na.value='grey75') + labs(x=NULL,y='Total mobile-DNA length per sample (kb)') + theme_classic(base_size=9) + theme(legend.position='none',axis.text.x=element_text(angle=30,hjust=1))
# B: element-level cargo proportion using independent element IDs
cargo <- read.delim(file.path(root,'cargo_records_final.tsv'),check.names=FALSE) %>% rename(cargo=category) %>% filter(cargo %in% c('AMR','STRESS','VFDB','defense')) %>% mutate(cargo=recode(cargo,STRESS='Stress',defense='Defence'))
# restrict to groups and exact element inventory; map cargo element IDs to inventory where possible
base <- inv %>% filter(BAPS %in% ord) %>% distinct(sample_id,BAPS,element_type,element_id)
# source records are deduplicated by sample/type/element/cargo; for records whose IDs differ, retain source element ID as unit
pos <- cargo %>% left_join(inv %>% distinct(sample_id,BAPS),by='sample_id') %>% filter(BAPS %in% ord) %>% distinct(sample_id,BAPS,element_type,element_id=cargo_feature,cargo)
den <- base %>% count(BAPS,element_type,name='n_elements')
num <- pos %>% count(BAPS,element_type,cargo,name='n_positive')
prev <- den %>% tidyr::expand_grid(cargo=c('AMR','Stress','VFDB','Defence')) %>% left_join(num,by=c('BAPS','element_type','cargo')) %>% mutate(n_positive=replace_na(n_positive,0),proportion=n_positive/n_elements,label=paste0(n_positive,'/',n_elements))
prev$BAPS <- factor(prev$BAPS,levels=ord)
pB <- ggplot(prev,aes(proportion,element_type,colour=cargo)) + geom_vline(xintercept=c(.25,.5,.75),colour='grey92') + geom_point(position=position_dodge(width=.58),size=2.5) + geom_text(aes(label=label),position=position_dodge(width=.58),hjust=-.15,size=2.3,show.legend=FALSE) + facet_wrap(~BAPS,nrow=1) + scale_x_continuous(labels=scales::percent,limits=c(0,1.12)) + scale_colour_brewer(palette='Dark2') + labs(x='Proportion of mobile elements carrying cargo',y=NULL,colour='Cargo') + theme_classic(base_size=9)
# C: defence systems per sample, by broad DefenseFinder type
D <- pos %>% filter(cargo=='Defence') %>% distinct(sample_id,BAPS,element_type,element_id,cargo) %>% filter(BAPS %in% ord) %>% count(sample_id,BAPS,name='defence_elements')
pC <- ggplot(D,aes(BAPS,defence_elements,fill=BAPS)) + geom_boxplot(outlier.shape=NA,width=.62,alpha=.25) + geom_jitter(width=.12,size=.7,alpha=.55) + scale_fill_manual(values=BAPS_PLOT_COLORS, na.value='grey75') + labs(x=NULL,y='Mobile-element defence systems per sample') + theme_classic(base_size=9) + theme(legend.position='none')
# C: original defence-system category composition by BAPS
Dall <- read.delim(file.path(root,'sample_mobilome_summary.tsv'),check.names=FALSE)
dc <- grep('^defense_',names(Dall),value=TRUE)
broad <- function(x){ z<-tolower(x); dplyr::case_when(grepl('cas',z)~'Cas/CRISPR', grepl('rm_',z)~'RM', grepl('brex',z)~'BREX', grepl('disarm',z)~'DISARM', grepl('cbass',z)~'CBASS', grepl('abi',z)~'Abi', grepl('retron',z)~'Retrons', grepl('maz|hip|hok|ta',z)~'Toxin–antitoxin', TRUE~'Other') }
origC <- Dall %>% select(sample_id,BAPS,all_of(dc)) %>% filter(BAPS %in% ord) %>% pivot_longer(all_of(dc),names_to='system',values_to='n') %>% mutate(system=broad(system),present=n>0) %>% group_by(BAPS,system) %>% summarise(prevalence=mean(present),.groups='drop')
origC$BAPS <- factor(origC$BAPS,levels=ord)
pC <- ggplot(origC,aes(system,prevalence,fill=BAPS)) + geom_col(position='dodge') + scale_y_continuous(labels=scales::percent) + scale_fill_manual(values=BAPS_PLOT_COLORS, na.value='grey75') + labs(x=NULL,y='Samples carrying defence-system category',fill='FastBAPS') + theme_classic(base_size=9) + theme(axis.text.x=element_text(angle=35,hjust=1))
# E: sample-level mobile-element defence burden
# D: bipartite network, edge weight = independent source element records
edge <- pos %>% count(element_type,cargo,name='weight')
nodes <- data.frame(name=unique(c(as.character(edge$element_type),as.character(edge$cargo))))
g <- graph_from_data_frame(edge,directed=FALSE,vertices=nodes); V(g)$kind <- ifelse(V(g)$name %in% unique(edge$element_type),'Mobile element','Cargo'); V(g)$type <- V(g)$kind=='Cargo'
pD <- ggraph(g,layout='bipartite') + geom_edge_link(aes(width=weight),colour='grey55',alpha=.7,show.legend=FALSE) + geom_node_point(aes(colour=kind),size=5) + geom_node_text(aes(label=name),size=2.8,repel=TRUE) + scale_colour_manual(values=c('Mobile element'='#1B9E77','Cargo'='#D95F02')) + theme_void() + labs(title='Mobile-element–cargo associations')
pE <- ggplot(D,aes(BAPS,defence_elements,fill=BAPS)) + geom_boxplot(outlier.shape=NA,width=.62,alpha=.25) + geom_jitter(width=.12,size=.7,alpha=.55) + scale_fill_manual(values=BAPS_PLOT_COLORS, na.value='grey75') + labs(x=NULL,y='Mobile-element defence systems per sample') + theme_classic(base_size=9) + theme(legend.position='none')
# F: relationship between mobile-DNA length and defence burden
length_by_sample <- len %>% group_by(sample_id,BAPS,element_type) %>% summarise(total_length_kb=sum(total_length_kb),.groups='drop')
defence_by_sample <- D %>% select(sample_id,BAPS,defence_elements)
rel <- length_by_sample %>% inner_join(defence_by_sample,by=c('sample_id','BAPS')) %>% filter(BAPS %in% ord, element_type != 'compositional_outlier')
cor_stats <- rel %>% group_by(element_type) %>% summarise(n=n(),spearman_rho= suppressWarnings(cor(total_length_kb,defence_elements,method='spearman',use='complete.obs')),p_value= suppressWarnings(cor.test(total_length_kb,defence_elements,method='spearman',exact=FALSE)$p.value),.groups='drop')
write.table(cor_stats,file.path(out,'Figure5F_length_defence_spearman.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
ann <- cor_stats %>% filter(element_type != 'compositional_outlier') %>% mutate(label=sprintf('R = %.2f',spearman_rho))
pF <- ggplot(rel,aes(total_length_kb,defence_elements,colour=BAPS)) + geom_point(alpha=.6,size=1.2) + geom_smooth(method='lm',se=TRUE,colour='grey35',fill='grey75',linewidth=.45) + facet_wrap(~element_type,scales='free_x',ncol=2) + geom_text(data=ann,aes(x=Inf,y=Inf,label=label),inherit.aes=FALSE,hjust=1.1,vjust=1.4,colour='black',size=3) + scale_colour_manual(values=BAPS_PLOT_COLORS, na.value='grey75') + labs(x='Total mobile-DNA length per sample (kb)',y='Mobile-element defence systems per sample',colour='FastBAPS') + theme_classic(base_size=9)

fig <- pA / (pB|pC) / (pD|pE) / pF + plot_annotation(tag_levels='A')
png(file.path(out,'Figure5_mobile_defence_cargo_final.png'),width=4200,height=4200,res=300); print(fig); dev.off(); ggsave(file.path(out,'Figure5_mobile_defence_cargo_final.pdf'),fig,width=14,height=14)
write.table(len,file.path(out,'Figure5A_sample_total_lengths.tsv'),sep='\t',row.names=FALSE,quote=FALSE); write.table(prev,file.path(out,'Figure5B_element_cargo_prevalence.tsv'),sep='\t',row.names=FALSE,quote=FALSE); write.table(edge,file.path(out,'Figure5D_network_edges.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
