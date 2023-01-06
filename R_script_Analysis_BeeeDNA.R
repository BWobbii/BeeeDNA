##### 2022 Sickel et al - eDNA wild bees

#load required packages
#data manipulation
library(tidyverse)
library(phyloseq) 

#statistics
library(rstatix)
library(PMCMRplus)
library(broom)
library(parsnip)
library(yardstick)
library(RVAideMemoire)


#visualisation
library(cowplot)
library(EnvStats)
library(ggsignif)


#data import
m.tib <- read_delim('2020-12_SampleFileCompl.txt', delim = '\t', col_names = T, col_select = 1:26)

#DNA quantity and quality

#keep experiment 1 samples only
#remove sample duplicates
#convert from character to vector

M.COI.yield <- m.tib %>% filter(Gene == unique(Gene[1]) & 
					Exp.Group == 'normal' &
					Replicate == 'A') %>%
			arrange(SampleID.DNA) %>%
			mutate(Group.Sub = factor(Group.Sub, levels = c('1-5', '6-10', '10+'))) %>%
			mutate_at(c('Exp.Group', 'SampleType'), as.factor)


#DNA yield, coi.short fragment
M.COI.yield %>% get_summary_stats(gDNA.mass.ng, type = "common") %>%
	as.data.frame()

#mean DNA yield & SD for fecal samples & swab samples separately
M.COI.yield %>% group_by(SampleType) %>%
	get_summary_stats(gDNA.mass.ng, type = "common") %>%
	as.data.frame()



#histogram, test for normal distribution
hist(M.COI.yield$gDNA.mass.ng)
shapiro.test(M.COI.yield$gDNA.mass.ng) #p < 0.05; not normally distributed


#Kruskal-Wallis test for sample type comparison
M.COI.yield %>% kruskal_test(gDNA.mass.ng ~ SampleType)
# <chr>        <int>     <dbl> <int> <dbl> <chr>         
# gDNA.mass.ng   100     0.495     1 0.482 Kruskal-Wallis	

#nest size
M.COI.yield %>% kruskal_test(gDNA.mass.ng ~ Group.Sub)

#post-hoc test, Kruskal-Nemenyi's All-Pairs Rank Comparison test
M.COI.yield %>% kwAllPairsNemenyiTest(gDNA.mass.ng ~ Group.Sub, data = .)

#kruskal test for sample types separately, ~ nest size
M.COI.yield %>% split(.$SampleType) %>%
	map(~kruskal_test(gDNA.mass.ng ~ Group.Sub, data = .))		

#add kwAllPairsNemenyiTest
M.COI.yield %>% split(.$SampleType) %>%
	map(~kwAllPairsNemenyiTest(gDNA.mass.ng ~ Group.Sub, data = .))		



#VISUALISATION, FIG 1

#custom function
yield_plot <- function(data, x, y) {
  ggplot(data, aes({{x}}, {{y}})) +
    geom_violin(alpha = .5, fill = 'gray')+
	geom_jitter(size = 2, pch = 21, fill = 'gray')+
	theme_bw() + theme(panel.grid.major=element_blank(),
		panel.grid.minor=element_blank(), axis.line=element_line(colour = 'black'))+
		 stat_n_text(geom="label", show.legend = FALSE)+
				geom_signif(comparisons = list(c("1-5", "10+")), annotations='*' ,
              map_signif_level=TRUE, y = 2010) +
		geom_signif(comparisons = list(c("1-5", "6-10")), annotations='*' ,
              map_signif_level=TRUE)
	
}


BP.yield1 <- yield_plot(M.COI.yield, Group.Sub, gDNA.mass.ng) +
	ylab('DNA yield [ng]') +
	theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

BP.yield1

BP.yieldS <- M.COI.yield %>% filter(SampleType == "swab") %>%
	yield_plot(.,Group.Sub,gDNA.mass.ng ) + 
	ylab('DNA yield [ng]') +
		theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

BP.yieldS

BP.yieldF <-  M.COI.yield %>% filter(SampleType == "feces") %>%
	yield_plot(.,Group.Sub,gDNA.mass.ng )+
	xlab('nest size') + ylab('DNA yield [ng]') +
	 scale_x_discrete(labels = c('S', 'M', 'L'))

BP.yieldF

plot_grid(BP.yield1, BP.yieldS, BP.yieldF, ncol = 1, labels = c('A', 'B', 'C'))
ggsave('YIELD.pdf', height = 10)
ggsave('YIELD.jpg', height = 10)



#DNA QUALITY ANALYSIS
#coi.short fragment, experiment 1
M.COI.sh <- filter(m.tib, Gene == unique(Gene)[1] & 
					Exp.Group == 'normal') %>%
			arrange(SampleID.DNA) 


M.sh.spl <- M.COI.sh %>%  split(.$Replicate)

#rearrange data to analyse PCR results between sample duplicates
PCR.sh <- M.sh.spl$A %>%
	rename(BandA = Band.1stPCR) %>%
	mutate(BandB = M.sh.spl$B$Band.1stPCR,
	DiffPCR = ifelse(BandA - BandB == 0, 1,0)) 

#custom functions
pcr_pos <- function(data, x){
	data %>%
		summarise(POS = sum({{x}}), TOT = length({{x}})) %>%
	mutate(PCT = POS/TOT, NEG = TOT-POS)
		}

pcr_agree <- function(data, x){
	data %>%
		summarise(AGREE = sum({{x}}),
				DISAGREE = length({{x}}[{{x}} ==0]),
				TOT = length({{x}})) %>%
		mutate(PCR.AGREE = AGREE/TOT)
		}

test_fisher <- function(data, var1, var2){
	data %>%	summarise(Y = {{var1}}, N = {{var2}}) %>%
			select(c(Y, N)) %>%
			as.data.frame() %>% fisher.test(alternative ="two.sided") 
			}


M.COI.sh %>% pcr_pos(Band.1stPCR)
PCR.sh %>% pcr_agree(DiffPCR)

#analysis by sample type
M.COI.sh %>% group_by(SampleType) %>% pcr_pos(Band.1stPCR)
PCR.sh %>% group_by(SampleType) %>% pcr_agree(DiffPCR)



M.COI.sh %>% group_by(SampleType) %>% pcr_pos(Band.1stPCR) %>%
		test_fisher(., POS, NEG) 

PCR.sh %>% group_by(SampleType) %>% pcr_agree(DiffPCR) %>%
		test_fisher(., AGREE, DISAGREE) 

#analysis by nest size, data for supplement
#both sample types
M.COI.sh %>% group_by(Group.Sub) %>% pcr_pos(Band.1stPCR)
PCR.sh %>% group_by(Group.Sub) %>% pcr_agree(DiffPCR)


M.COI.sh %>% group_by(Group.Sub) %>% pcr_pos(Band.1stPCR) %>%
		test_fisher(., POS, NEG) #not significant

PCR.sh %>% group_by(Group.Sub) %>% pcr_agree(DiffPCR) %>%
		test_fisher(., AGREE, DISAGREE) 


PCR.sh %>% group_by(Group.Sub) %>% pcr_agree(DiffPCR) %>%
	column_to_rownames(var = "Group.Sub") %>%
	select(AGREE, DISAGREE) %>% 	t() %>% 
	fisher.multcomp




#nest size & sample types separately
M.COI.sh %>% filter(SampleType == 'feces') %>% group_by(Group.Sub) %>% pcr_pos(Band.1stPCR)
PCR.sh %>% filter(SampleType == 'feces') %>% group_by(Group.Sub) %>% pcr_agree(DiffPCR)

M.COI.sh %>% filter(SampleType == 'feces') %>% group_by(Group.Sub) %>% pcr_pos(Band.1stPCR) %>%
		test_fisher(., POS, NEG) 

PCR.sh %>% filter(SampleType == 'feces') %>% group_by(Group.Sub) %>% pcr_agree(DiffPCR) %>%
		test_fisher(., AGREE, DISAGREE) 



PCR.sh %>% filter(SampleType == 'feces') %>% group_by(Group.Sub) %>% pcr_agree(DiffPCR) %>%
	column_to_rownames(var = "Group.Sub") %>%
	select(AGREE, DISAGREE) %>% 	t() %>% 
	fisher.multcomp




M.COI.sh %>% filter(SampleType == 'swab') %>% group_by(Group.Sub) %>% pcr_pos(Band.1stPCR)
PCR.sh %>% filter(SampleType == 'swab') %>% group_by(Group.Sub) %>% pcr_agree(DiffPCR)

M.COI.sh %>% filter(SampleType == 'swab') %>% group_by(Group.Sub) %>% pcr_pos(Band.1stPCR) %>%
		test_fisher(., POS, NEG)

PCR.sh %>% filter(SampleType == 'swab') %>% group_by(Group.Sub) %>% pcr_agree(DiffPCR) %>%
		test_fisher(., AGREE, DISAGREE) 

PCR.sh %>% filter(SampleType == 'swab') %>% group_by(Group.Sub) %>% pcr_agree(DiffPCR) %>%
	column_to_rownames(var = "Group.Sub") %>%
	select(AGREE, DISAGREE) %>% 	t() %>% 
	fisher.multcomp



#coi.long fragment
M.COI.lg <- m.tib %>% filter(Gene == unique(Gene)[2] & 
					Exp.Group == 'normal') %>%
			arrange(SampleID.DNA) 

M.lg.spl <- M.COI.lg %>% split(.$Replicate)


#rearrange data to analyse PCR results between sample duplicates
PCR.lg <- M.lg.spl$C %>%
	rename(BandC = Band.1stPCR) %>%
	mutate(BandD = M.lg.spl$D$Band.1stPCR,
	DiffPCR = ifelse(BandC - BandD == 0, 1,0)) 



#PCR positives, all samples
M.COI.lg %>% pcr_pos(Band.1stPCR) 
PCR.lg %>% pcr_agree(DiffPCR)


#analysis by nest size
M.COI.lg %>% group_by(Group.Sub) %>% pcr_pos(Band.1stPCR)
PCR.lg %>% group_by(Group.Sub) %>% pcr_agree(DiffPCR)

M.COI.lg %>% group_by(Group.Sub) %>% pcr_pos(Band.1stPCR) %>%
		test_fisher(., POS, NEG) 
PCR.lg %>% group_by(Group.Sub) %>% pcr_agree(DiffPCR) %>%
		test_fisher(., AGREE, DISAGREE) 




#### SEQUENCING DATA, COI.SHORT FRAGMENT
m.seq <- m.tib %>% filter(Gene == unique(m.tib$Gene)[1]) %>%
		mutate(NEST.ID = paste(Trap.ID, Nest.No, sep = '.'),
			Exp.Group = factor(ifelse(Exp.Group == 'mixed', 'mixed', 
							ifelse(Exp.Group == 'para', 'para', 
								ifelse(Exp.Group == 'normal' & SampleType == 'swab', 'swab','feces'))),
						levels = c('feces', 'swab', 'para', 'mixed')),

				Group.Sub = factor(ifelse(Exp.Group %in% c('swab', 'feces'),Group.Sub, 
							ifelse(No.Cells <= 5, '1-5',
								ifelse(No.Cells <= 10, '6-10', '10+'))),
						levels = c('1-5', '6-10', '10+'))) %>%
		as.data.frame()
 

rownames(m.seq) <- m.seq[,1]

#import ASV table
ASV.sh <- import_biom('asvtab.biom')

#import taxonomy table
TAX.V <- tax_table(as.matrix(read.table("taxonomy.vsearch", header=T,row.names=1,fill=T, sep = ",")))

#merge into phyloseq object 
COI.sh <- phyloseq(otu_table(ASV.sh, taxa_are_rows=TRUE), tax_table(TAX.V), sample_data(m.seq))

#subset to Metazoa
(COI.meta <- subset_taxa(COI.sh, kingdom == "k:Metazoa")) #274 taxa

#analyse laboratory controls
(controls <- subset_samples(COI.meta, SampleType %in% c("PCRControl", "ExtrControl")))
(controls		= prune_taxa(taxa_sums(controls)>75, controls))

#remove contaminating taxa from dataset
(cont.ASVs <- rownames(tax_table(controls)))
(COI.meta <- subset_taxa(COI.meta, !taxa_names(COI.meta) %in% cont.ASVs))
#257 taxa

#subset samples plate-wise
subset_samples(COI.meta, SampleType == "PosControl") %>%
	subset_taxa(species == "s:Apis_mellifera") %>%
	otu_table()



#subset to samples only, remove Apis
COI.meta <- subset_samples(COI.meta, SampleType %in% c("swab", "feces")) %>%
	subset_taxa(genus != "g:Apis")
#256 taxa

(COI.sp <- tax_glom(COI.meta, "species")) #55 taxa
(tax_glom(COI.meta, "genus")) #49 taxa
(tax_glom(COI.meta, "family")) #31 taxa

#subset to Arthropods
#then to Hymenoptera + Diptera only
#then to Hymenoptera only
(sh.Arthro <- subset_taxa(COI.sp, phylum == "p:Arthropoda")) #55 taxa
(sh.HyDip <- subset_taxa(COI.sp, order %in% c("o:Hymenoptera", "o:Diptera"))) #17 taxa
(sh.Hy.all <- subset_taxa(COI.sp, order == "o:Hymenoptera")) #13 taxa

#number of Hymenoptera detected
sum(as.numeric(colSums(otu_table(sh.Hy.all)) > 0))

#remove unclassified species
(sh.Hy <- subset_taxa(sh.Hy.all, species != "s:_")) #12 taxa

COI.sp %>% estimate_richness(measures = 'Observed') %>%
			get_summary_stats(type='common') %>%
			as.data.frame()


sh.HyDip %>% estimate_richness(measures = 'Observed') %>%
			get_summary_stats(type='common') %>%
			as.data.frame()


#counting Hymenoptera detections by converting ASV-data from sh.Hy-object
#counting mixed species dections by converting richness data from sh.HyDip-object
as.numeric(colSums(otu_table(sh.Hy)) > 0) %>% sum

HY.tib <- as_tibble(data.frame(sample_data(sh.HyDip)[,c(1:15,25)])) %>%
	mutate(Hy.det = as.numeric(colSums(otu_table(sh.Hy)) > 0),
		N.spec = as.numeric(estimate_richness(sh.HyDip, measures="Observed")[,1] > 1),
		SampleType = factor(SampleType, , levels = c('feces', 'swab'))) %>%
		mutate_at(c('Replicate', 'Group.Sub', 'Exp.Group'), factor)

HY.DET <- HY.tib %>%
	group_by(Exp.Group, Group.Sub, SampleType) %>%
		summarise(n.Hy = sum(Hy.det == 1),
			total = n(), n.NoHy = (total- n.Hy)) %>%
	mutate(pct = round((n.Hy/total *100), digits =1))



#VISUALISATION, FIG 3

#custom function
hydet_plot <- function(data, x, y, label) {
  ggplot(data, aes({{x}}, {{y}})) +
    geom_col(position = position_dodge(), colour = 'black') +
	ylab('Samples with Hymenoptera \n at species level [%]')+
	theme_bw() + theme(panel.grid.major=element_blank(),
		panel.grid.minor=element_blank(), axis.line=element_line(colour = 'black'))+
	geom_text(aes({{x}}, {{y}}, label = {{label}}), position = position_dodge(width = 1), vjust = -0.5)+
	ylim(-0.5,100) +
	theme(legend.position = "top") +
	geom_hline(yintercept = 50) 
		}

#hym det by sample type
P.Hydet1 <- HY.DET %>% group_by(SampleType) %>% 
		summarise(n.HySp = sum(n.Hy), tot = sum(total)) %>%
		mutate(pct = round((n.HySp/tot*100), digits = 1)) %>%
		hydet_plot(., SampleType, pct, n.HySp) +
		xlab('Sample type')+
		scale_x_discrete(labels = c('feces', 'swab (exp. 1&2)'))+
		geom_signif(comparisons = list(c("feces", "swab")), annotations='*' ,
              map_signif_level=TRUE, y = 75) +
		geom_text(aes(x = SampleType, y = 95, label = c('N = 100', 'N = 210')), position = position_dodge(width = 1), vjust = -0.5, size = 4)



#hym det by experimental group, swabs only
P.Hydet2 <- HY.DET %>% filter(SampleType == "swab") %>% group_by(Exp.Group) %>%
	summarise(n.HySp = sum(n.Hy), tot = sum(total)) %>%
	mutate(pct = round((n.HySp/tot*100), digits = 1)) %>%
		hydet_plot(., Exp.Group, pct, n.HySp) +
		xlab('experimental group') +
		geom_text(aes(x = Exp.Group, y = 95, label = c('N = 100', 'N = 74', 'N = 36')), position = position_dodge(width = 1), vjust = -0.5, size = 4)+
		scale_fill_discrete(name = "Experimental group")



plot_grid(P.Hydet1, P.Hydet2, ncol =1, labels = c('A', 'B'))
ggsave('HY.DET.OCT.pdf', height = 7.5)
ggsave('HY.DET.OCT.jpg', height = 7.5)




#comparison between sample types & nest sizes
#BY SAMPLE TYPE / EXPERIMENTAL GROUP
HY.DET %>% group_by(SampleType) %>% 
		summarise(n.HySp = sum(n.Hy), tot = sum(total)) %>%
		mutate(pct = round((n.HySp/tot*100), digits = 1))

HY.DET %>% filter(SampleType == "swab") %>% group_by(Exp.Group) %>%
	summarise(n.HySp = sum(n.Hy), tot = sum(total)) %>%
	mutate(pct = round((n.HySp/tot*100), digits = 1))



test_fisher <- function(data, group, var1, var2){
	data %>% group_by({{group}}) %>% 
			summarise(Y = sum({{var1}}), N = sum({{var2}})) %>%
			select(c(Y, N)) %>%
			as.data.frame() %>% fisher.test(alternative ="two.sided") 
			}

#sample type, incl. para&mixed
HY.DET %>% test_fisher(., SampleType, n.Hy, n.NoHy)


#swabs only, by exp. group
HY.DET %>% filter(SampleType == 'swab') %>%
		test_fisher(., Exp.Group, n.Hy, n.NoHy)


#incl. unclassified Hymenoptera
HY.tib.unclass <- as_tibble(data.frame(sample_data(sh.HyDip)[,c(1:15,25)])) %>%
	mutate(Hy.det = as.numeric(colSums(otu_table(sh.Hy.all)) > 0),
		N.spec = as.numeric(estimate_richness(sh.HyDip, measures="Observed")[,1] > 1),
		SampleType = factor(SampleType, , levels = c('feces', 'swab'))) %>%
		mutate_at(c('Replicate', 'Group.Sub', 'Exp.Group'), factor)


HY.DET.unclass <- HY.tib.unclass %>%
	group_by(Exp.Group, Group.Sub, SampleType) %>%
		summarise(n.Hy = sum(Hy.det == 1),
			total = n(), n.NoHy = (total- n.Hy)) %>%
	mutate(pct = round((n.Hy/total *100), digits =1))

HY.DET.unclass %>% group_by(SampleType) %>% 
		summarise(n.HySp = sum(n.Hy), tot = sum(total)) %>%
		mutate(pct = round((n.HySp/tot*100), digits = 1))

HY.DET.unclass %>% filter(SampleType == "swab") %>% group_by(Exp.Group) %>%
	summarise(n.HySp = sum(n.Hy), tot = sum(total)) %>%
	mutate(pct = round((n.HySp/tot*100), digits = 1))


#sample type, incl. para&mixed
HY.DET.unclass %>% test_fisher(., SampleType, n.Hy, n.NoHy)


#swabs only, by exp. group
HY.DET.unclass %>% filter(SampleType == 'swab') %>%
		test_fisher(., Exp.Group, n.Hy, n.NoHy)




#### BY NEST SIZE
HY.DET %>% 	test_fisher(., Group.Sub, n.Hy, n.NoHy)


#REPEAT FOR SAMPLE TYPE SEPARATELY (NOT SHOWN)
HY.DET %>% filter(SampleType == 'feces') %>%
			test_fisher(., Group.Sub, n.Hy, n.NoHy)

HY.DET %>% filter(SampleType == 'swab') %>%
			test_fisher(., Group.Sub, n.Hy, n.NoHy)


#####################################

#GLM ANALYSIS
#function to reshape the data

reshape_det2 <- function(data, replicate1, replicate2){
	data %>%
	filter(Replicate == replicate1) %>%
	rename(Hy1 = Hy.det) %>%
	mutate(Hy2 = data$Hy.det[data$Replicate == replicate2],
		DiffHy = ifelse(Hy1 + Hy2 == 2, 1, 0)) %>%
	mutate_at(c('SampleType', 'Group.Sub'), as.numeric)
	}

reshape_det1 <- function(data, replicate1, replicate2){
	data %>%
	filter(Replicate == replicate1) %>%
	rename(Hy1 = Hy.det) %>%
	mutate(Hy2 = data$Hy.det[data$Replicate == replicate2],
		DiffHy = ifelse(Hy1 + Hy2 > 0, 1, 0)) %>%
	mutate_at(c('SampleType', 'Group.Sub'), as.numeric)
	}



GLM_details <- logistic_reg() %>% set_engine("glm", family = binomial(link = 'logit')) %>%
			set_mode("classification") %>%
				 translate()



GLM_analysis1 <- function(model, null.model) {
		GLM_res1 <- glance(model) %>%								#goodness of fit, pseudo R sq (McFadden)
		  mutate(pchisq = pchisq(deviance, df.residual, lower.tail = FALSE),
				pseudo.rsq = 1 - deviance/null.deviance) %>%
  				select(pchisq, pseudo.rsq) 

		LRT <- anova(								#LR test on the model
  			extract_fit_engine(null.model),
			extract_fit_engine(model),
 			 test = "LRT") 
		LRT <- LRT$Pr[2]
		
		Wald <- Anova(							#TypeII Wald Test					
			extract_fit_engine(model), type = "II", test = "Wald") %>% 
				tidy()
		
		GLM_res1 <- bind_cols(GLM_res1, LRT = LRT,Wald)
		return(GLM_res1)
			}


#
GLM_analysis2 <- function(model, data, var){
		pred <- model %>% predict(new_data = data) %>% mutate(truth = var)					#predict
		conf_mat <- conf_mat(pred, truth = truth, estimate = .pred_class) %>% 
		tidy %>% mutate(name = c('TN', 'FP', 'FN', 'TP'))%>% list
		acc <- accuracy(pred, truth = truth, estimate = .pred_class) #yes
		prec <- precision(pred, truth = truth, estimate = .pred_class) #yes
		recall <- 	recall(pred, truth = truth, estimate = .pred_class)
		F1 <- 	f_meas(pred, truth = truth, estimate = .pred_class)
		GLM_res <- tibble(acc = acc$.estimate, 
					prec = prec$.estimate, 
					recall = recall$.estimate, 
					F1 = F1$.estimate, 
					conf_mat = conf_mat)
		return(GLM_res)
			}



#GLM No.1: HYM DETECTION IN BOTH DUPLICATES
HY.DET2 <- HY.tib %>% filter(Exp.Group %in% c('feces', 'swab')) %>% reshape_det2(., 'A', 'B') %>% 
			select(DiffHy, gDNA.mass.ng, SampleType, Group.Sub)

HY.DET2 %>%	summarise(count = sum(DiffHy),
			pct = sum(DiffHy)/length(DiffHy))	

HY.DET2 <- HY.DET2 %>% mutate_at('DiffHy', factor)


HY.DET2_FIT <-  GLM_details %>% fit(DiffHy ~ gDNA.mass.ng * SampleType + gDNA.mass.ng * Group.Sub, data = HY.DET2)
HY.DET2_NULL <-  GLM_details %>% fit(DiffHy ~ 1, data = HY.DET2)

RES.HY.DET2_1 <- GLM_analysis1(model = HY.DET2_FIT, null.model = HY.DET2_NULL )
RES.HY.DET2_2 <- GLM_analysis2(model = HY.DET2_FIT, data = HY.DET2, var = HY.DET2$DiffHy)

RES.HY.DET2_1
RES.HY.DET2_2
RES.HY.DET2_2$conf_mat

#GLM No.2: HYM DETECTION IN AT LEAST ONE REPLICATE

HY.DET1 <-HY.tib %>% filter(Exp.Group %in% c('feces', 'swab')) %>% reshape_det1(., 'A', 'B') %>% 
			select(DiffHy, gDNA.mass.ng, SampleType, Group.Sub)

HY.DET1 %>%	summarise(count = sum(DiffHy),
			pct = sum(DiffHy)/length(DiffHy))	

HY.DET1 <- HY.DET1 %>% mutate_at('DiffHy', factor)


HY.DET1_FIT <-  GLM_details %>% fit(DiffHy ~ gDNA.mass.ng * SampleType + gDNA.mass.ng * Group.Sub, data = HY.DET1)
HY.DET1_NULL <-  GLM_details %>% fit(DiffHy ~ 1, data = HY.DET1)

RES.HY.DET1_1 <- GLM_analysis1(model = HY.DET1_FIT, null.model = HY.DET1_NULL )
RES.HY.DET1_2 <- GLM_analysis2(model = HY.DET1_FIT, data = HY.DET1, var = HY.DET1$DiffHy)

RES.HY.DET1_1
RES.HY.DET1_2
RES.HY.DET1_2$conf_mat



#GLM No.3: HYM DETECTION ~ PCR SUCCESS
#convert factors to numeric
HY.DETpcr <- HY.tib %>% mutate_at(c('SampleType', 'Group.Sub'), as.numeric)

HY.DETpcr %>%	summarise(count = sum(Hy.det),
			pct = sum(Hy.det)/length(Hy.det))	


HY.DETpcr <- HY.DETpcr %>% mutate_at('Hy.det', factor)

HY.DETpcr_FIT <-  GLM_details %>% fit(Hy.det ~ Band.1stPCR * SampleType + Band.1stPCR * Group.Sub, data = HY.DETpcr)
HY.DETpcr_NULL <-  GLM_details %>% fit(Hy.det ~ 1, data = HY.DETpcr)

RES.HY.DETpcr_1 <- GLM_analysis1(model = HY.DETpcr_FIT, null.model = HY.DETpcr_NULL )
RES.HY.DETpcr_2 <- GLM_analysis2(model = HY.DETpcr_FIT, data = HY.DETpcr, var = HY.DETpcr$Hy.det)


RES.HY.DETpcr_1
RES.HY.DETpcr_2
RES.HY.DETpcr_2$conf_mat




### MIXED SPECIES DETECTION, EXPERIMENT 2
HY.MIX <- filter(HY.tib, Exp.Group %in% c('para', 'mixed'))

HY.MIX %>%	summarise(count = sum(N.spec),
			pct = count/length(N.spec))	



HY.MIX %>% group_by(Exp.Group) %>%
	summarise(MIX = sum(N.spec), TOT = length(N.spec)) %>%
	mutate(PCT = MIX/TOT)

#for GLM analysis

reshape_mix2 <- function(data, replicate1, replicate2){
	data %>%
	filter(Replicate == replicate1) %>%
	rename(Mix1 = N.spec) %>%
	mutate(Mix2 = data$N.spec[data$Replicate == replicate2],
		DiffMix = ifelse(Mix1 + Mix2 == 2, 1, 0)) %>%
	mutate_at(c('SampleType', 'Group.Sub'), as.numeric)
	}

reshape_mix1 <- function(data, replicate1, replicate2){
	data %>%
	filter(Replicate == replicate1) %>%
	rename(Mix1 = N.spec) %>%
	mutate(Mix2 = data$N.spec[data$Replicate == replicate2],
		DiffMix = ifelse(Mix1 + Mix2 > 0, 1, 0)) %>%
	mutate_at(c('SampleType', 'Group.Sub'), as.numeric)
	}


#GLM No.4: SPECIES MIX IN BOTH DUPLICATES
HY.MIX2 <- reshape_mix2(HY.MIX, 'A', 'B') %>% 
			select(DiffMix, gDNA.mass.ng, SampleType, Group.Sub)

HY.MIX2 %>%	summarise(count = sum(DiffMix),
			pct = count/length(DiffMix))	

HY.MIX2 <- HY.MIX2 %>% mutate_at('DiffMix', factor)


HY.MIX2_FIT <-  GLM_details %>% fit(DiffMix ~ gDNA.mass.ng, data = HY.MIX2)
HY.MIX2_NULL <-  GLM_details %>% fit(DiffMix ~ 1, data = HY.MIX2)

RES.HY.MIX2_1 <- GLM_analysis1(model = HY.MIX2_FIT, null.model = HY.MIX2_NULL )
RES.HY.MIX2_2 <-GLM_analysis2(model = HY.MIX2_FIT, data = HY.MIX2, var = HY.MIX2$DiffMix)


RES.HY.MIX2_1
RES.HY.MIX2_2
RES.HY.MIX2_2$conf_mat



#GLM No.5: SPECIES MIX IN AT LEAST ONE REPLICATE
HY.MIX1 <- reshape_mix1(HY.MIX, 'A', 'B') %>% 
			select(DiffMix, gDNA.mass.ng, SampleType, Group.Sub)

HY.MIX1 %>%	summarise(count = sum(DiffMix),
			pct = count/length(DiffMix))	

HY.MIX1 <- HY.MIX1 %>% mutate_at('DiffMix', factor)


HY.MIX1_FIT <-  GLM_details %>% fit(DiffMix ~ gDNA.mass.ng, data = HY.MIX1)
HY.MIX1_NULL <-  GLM_details %>% fit(DiffMix ~ 1, data = HY.MIX1)

RES.HY.MIX1_1 <- GLM_analysis1(model = HY.MIX1_FIT, null.model = HY.MIX1_NULL )
RES.HY.MIX1_2 <- GLM_analysis2(model = HY.MIX1_FIT, data = HY.MIX1, var = HY.MIX1$DiffMix)


RES.HY.MIX1_1
RES.HY.MIX1_2
RES.HY.MIX1_2$conf_mat



#GLM No.6: SPECIES MIX ~ PCR SUCCESS
#convert factors to numeric
HY.MIXpcr <- HY.MIX %>% mutate_at(c('SampleType', 'Group.Sub'), as.numeric)

HY.MIXpcr %>%	summarise(count = sum(N.spec),
			pct = count/length(N.spec))	


HY.MIXpcr <- HY.MIXpcr %>% mutate_at('N.spec', factor)

HY.MIXpcr_FIT <-  GLM_details %>% fit(N.spec ~ Band.1stPCR, data = HY.MIXpcr)
HY.MIXpcr_NULL <-  GLM_details %>% fit(N.spec ~ 1, data = HY.MIXpcr)

RES.HY.MIXpcr_1 <- GLM_analysis1(model = HY.MIXpcr_FIT, null.model = HY.MIXpcr_NULL )
RES.HY.MIXpcr_2 <- GLM_analysis2(model = HY.MIXpcr_FIT, data = HY.MIXpcr, var = HY.MIXpcr$N.spec)


RES.HY.MIXpcr_1
RES.HY.MIXpcr_2
RES.HY.MIXpcr_2$conf_mat




### eDNA data validation with morphological data 
#import morpholigal/visual data

MORPHO <-  read_delim('2022-09_RESULTS_MORPHO.txt', delim = '\t', col_names = T)



#number of nests, for which morphological analysis was available
MORPHO$species1_morpho %>% na.omit %>% length #57

#use eDNA data of all detected Arthropods
DNA.arthro <- data.frame(otu_table(sh.Arthro))
DNA.tax <- data.frame(tax_table(sh.Arthro))

#number of samples without any Arthropod detection
(no.DNA <- length(as.numeric(colSums(DNA.arthro) > 0)) - sum(as.numeric(colSums(DNA.arthro) > 0)))
#7

#overall comparion between morphological data and eDNA
sp.DNA <- get_taxa_unique(sh.Arthro, "species")
sp.DNA <- subset(sp.DNA, sp.DNA != "s:_")
sp.MORPHO <- MORPHO %>% pull(species1_morpho, species2_morpho) %>% unique %>% na.omit

(sp.T1 <- sp.MORPHO[sp.MORPHO %in% sp.DNA == "TRUE"])
#in agreement: 3 wild bees (O. bicornis, H. communis, H. truncorum), 
#2 solitary wasps (Trypoxylon sp.), 
#2 parasitods (C. indagaor, S. decemguttata)


(sp.F1 <- sp.MORPHO[sp.MORPHO %in% sp.DNA == "FALSE"])
#not detected via eDNA: H. leucomelana, T. cyanea

(sp.T2 <- sp.DNA[sp.DNA %in% sp.MORPHO == "TRUE"])
#as above
(sp.F2 <- sp.DNA[sp.DNA %in% sp.MORPHO == "FALSE" ])
length(sp.F2)
#detected via eDNA, but not morphology - 46 taxa

sp.DNAonly <- subset_taxa(sh.Arthro, tax_table(sh.Arthro)[,7] %in% sp.F2)
get_taxa_unique(sp.DNAonly, 'order') 
write.table(tax_table(sp.DNAonly), 'taxa-DNAonly.txt')





DNA.taxa <- DNA.arthro %>% summarise(across(dplyr::everything(), 
			list(rowname = ~list(which(. > 0))))) %>%

	pivot_longer(
		cols = contains("_"),
		names_to = "key",
		values_to = "val",
		values_transform = list(val = as.character))	%>%


	separate(key, c("column", "name"), sep = "_") %>%

	pivot_wider(
		names_from = name,
		values_from = val) %>%

	mutate(rowname = str_replace_all(rowname, 
							c("c\\(" = "", 
								"\\)" = "",
								"integer\\(0" = "NA")),
		column = paste(column, get_variable(sh.Arthro, "NEST.ID"), sep = '-'),
		NEST.ID = get_variable(sh.Arthro, "NEST.ID"),
		Sample.ID = get_variable(sh.Arthro, "SampleID.ASV"),
		Sample.Type = get_variable(sh.Arthro, "SampleType"),
		Exp = get_variable(sh.Arthro, "Exp.Group"),
		Nest.Size = get_variable(sh.Arthro, "Group.Sub"))


DNA.taxa$zOTU <- lapply(DNA.taxa$rowname, function(x) as.numeric(strsplit(x, ",")[[1]])) %>%
		 rapply(., how = "replace", function(x) rownames(DNA.arthro)[x])


#custom function
replace_zOTU <- function(data.col, level) {
		rapply(data.col, how = "replace", function(x) DNA.tax[x,level]) 
		}


DNA.taxa$order_DNA <- replace_zOTU(DNA.taxa$zOTU, level = 'order')
DNA.taxa$family_DNA <- replace_zOTU(DNA.taxa$zOTU, level = 'family') %>%
	sapply( function(x) subset(x, x != "f:_")) 
DNA.taxa$genus_DNA <- replace_zOTU(DNA.taxa$zOTU, level = 'genus') %>%
	sapply( function(x) subset(x, x != "g:_")) 
DNA.taxa$species_DNA <- replace_zOTU(DNA.taxa$zOTU, level = 'species')  %>%
	sapply( function(x) subset(x, x != "s:_")) 



RES.COMP = NULL

for (i in 1:nrow(MORPHO)) {
DNA.taxa %>% 
	filter(NEST.ID == MORPHO$ID_nest[i]) %>% 
  mutate(DNA_morph_sp = map(species_DNA, str_detect, pattern = paste(c(MORPHO$species1_morpho[i], MORPHO$species2_morpho[i]), collapse = "|")),
         DNA_morph_sp = map_lgl(DNA_morph_sp, any),
		DNA_vis_sp = map(species_DNA, str_detect, pattern = paste(c(MORPHO$species1_visual[i], MORPHO$species2_visual[i]), collapse = "|")),
         DNA_vis_sp = map_lgl(DNA_vis_sp, any),

	DNA_morph_gen = map(genus_DNA, str_detect, pattern = paste(c(MORPHO$genus1_morpho[i], MORPHO$genus2_morpho[i]), collapse = "|")),
         DNA_morph_gen = map_lgl(DNA_morph_gen, any),
		
		DNA_vis_gen = map(genus_DNA, str_detect, pattern = paste(c(MORPHO$genus1_visual[i], MORPHO$genus2_visual[i]), collapse = "|")),
         DNA_vis_gen = map_lgl(DNA_vis_gen, any),
		

		DNA_morph_fam = map(family_DNA, str_detect, pattern = paste(c(MORPHO$family1_morpho[i], MORPHO$family2_morpho[i]), collapse = "|")),
         	DNA_morph_fam = map_lgl(DNA_morph_fam, any),
		
		DNA_vis_fam = map(family_DNA, str_detect, pattern = paste(c(MORPHO$family1_visual[i], MORPHO$family2_visual[i]), collapse = "|")),
         	DNA_vis_fam = map_lgl(DNA_vis_fam, any)		,


		DNA_morph_ord = map(order_DNA, str_detect, pattern = paste(c(MORPHO$order1_morpho[i], MORPHO$order2_morpho[i]), collapse = "|")),
         	DNA_morph_ord = map_lgl(DNA_morph_ord, any),
		
		DNA_vis_ord = map(order_DNA, str_detect, pattern = paste(c(MORPHO$order1_visual[i], MORPHO$order2_visual[i]), collapse = "|")),
         	DNA_vis_ord = map_lgl(DNA_vis_ord, any)) %>%

		select(3:7, 13:20) %>%
		mutate_if(is.logical, as.numeric) %>%
		bind_rows(., RES.COMP) -> RES.COMP
	}


#### Screening for
#Ichneumonoidea
#Chalicidoidea
#Spheciformes

MORPHO %>% filter(family2_morpho == "f:Ichneumonoidea") -> MORPHO.ichneu
MORPHO %>% filter(family2_visual == "f:Chalicidoidea") -> MORPHO.chali
MORPHO %>% filter(family1_visual == "f:Spheciformes") -> MORPHO.spheci
MORPHO %>% filter(family2_visual == "f:Spheciformes") %>% bind_rows(., MORPHO.spheci) -> MORPHO.spheci


TAX.ICHNEU <- list(family = "f:Ichneumonidae",
		genus = "g:Hoplocryptus|g:Poemenia|g:Stenarella|g:Perithous|g:Clistopyga",
		species = "s:Peomenia_collaris|s:Stenarella_domator|s:Perithous_scurra|s:Perithous_septemcinctorius|s:Perithous_divinator|s:Perithous_medicator|s:Clistopyga_incitator")

TAX.CHALI <- list(family = "f:Torymidae|f:Leucospidae|f:Eulophidae|f:Encyrtidae",
			genus = "g:Monodontomerus",
			species = "s:Monodontomerus_obscurus")


TAX.SPHECI <- list(family = "f:Crabronidae|f:Sphecidae",
			genus = "g:Passaloecus|g:Pemphredon|g:Psenulus|g:Trypoxylon",
			species = "s:Passaloecus_corniger|s:Passaloecus_gracilis|s:Passaloecus_insignis|s:Passaloecus_monilicornis|s:Passaloecus_singularis|s:Pemphredon_lethifer|s:Psenulus_chevrieri|s:Psenulus_concolor|s:Psenulus_pallipes|s:Trypoxylon_attenuatum|s:Trypoxylon_clavicerum|s:Trypoxylon_figulus|s:Trypoxylon_minus")


compare_special_tax <- function(data, morpho, taxa){
	
		data %>% filter(NEST.ID %in% morpho$ID_nest) %>%
				mutate(DNA_special_sp = map(species_DNA, str_detect, pattern = taxa$species),
					DNA_special_sp = map_lgl(DNA_special_sp, any),

					DNA_special_gen = map(genus_DNA, str_detect, pattern = taxa$genus),
					DNA_special_gen = map_lgl(DNA_special_gen, any),
			
				DNA_special_fam = map(family_DNA, str_detect, pattern = taxa$family),
					DNA_special_fam = map_lgl(DNA_special_fam, any)) %>%
				mutate_if(is.logical, as.numeric)
				}




RES.special <- compare_special_tax(DNA.taxa, MORPHO.ichneu,  TAX.ICHNEU)
RES.special <- compare_special_tax(DNA.taxa, MORPHO.chali,  TAX.CHALI) %>% bind_rows(., RES.special)
RES.special <- compare_special_tax(DNA.taxa, MORPHO.spheci,  TAX.SPHECI) %>% bind_rows(., RES.special)


RES.COMP <- RES.COMP %>% mutate_if(is.numeric,  ~replace(., is.na(.), 0))  %>% 
				mutate(comp.sp = ifelse(DNA_morph_sp + DNA_vis_sp > 0, 1,0),
				comp.gen = ifelse(DNA_morph_gen + DNA_vis_gen > 0, 1,0),
				comp.fam = ifelse(DNA_morph_fam + DNA_vis_fam > 0, 1,0),
				comp.ord = ifelse(DNA_morph_ord + DNA_vis_ord > 0, 1,0))


RES.add <- RES.special %>% filter(DNA_special_sp > 0 | DNA_special_gen > 0 | DNA_special_fam > 0)
RES.COMP$comp.sp[RES.COMP$Sample.ID %in% RES.add$Sample.ID] <- 1
RES.COMP$comp.gen[RES.COMP$Sample.ID %in% RES.add$Sample.ID] <- 1
RES.COMP$comp.fam[RES.COMP$Sample.ID %in% RES.add$Sample.ID] <- 1


RES.COMP2 <- tibble(Sample.ID = rep(RES.COMP$Sample.ID, 4),
			NEST.ID = rep(RES.COMP$NEST.ID, 4),
			Sample.Type = rep(RES.COMP$Sample.Type, 4),
			EXP = rep(RES.COMP$Exp, 4),
			NEST.Size = rep(RES.COMP$Nest.Size, 4),
			tax.level = factor(c(rep('species', 310),rep('genus', 310), rep('family', 310), rep('order', 310)),
							levels = c('order', 'family', 'genus', 'species')),
			ID.cor = c(RES.COMP$comp.sp,RES.COMP$comp.gen,RES.COMP$comp.fam, RES.COMP$comp.ord ))

#determine actual number of comparisons
#depending on entries in MORPHO-dataset
#and substracting 7 (= number of samples without an entry in DNA-dataset)

NEST_wSp <- MORPHO %>% 
	mutate(tmp = paste0(species1_morpho, species2_morpho, species1_visual, species2_visual),
		tmp2 = str_detect(tmp, 's:')) %>% filter(tmp2 == TRUE) %>% pull(ID_nest) 

TOT.SP <- (RES.COMP$NEST.ID %in% NEST_wSp %>% sum) -7


NEST_wGen <- MORPHO %>% 
mutate( tmp = paste0(genus1_morpho, genus2_morpho, genus1_visual, genus2_visual),
		tmp2 = str_detect(tmp, 'g:')) %>% filter(tmp2 == TRUE) %>% pull(ID_nest) 

TOT.GEN <- (RES.COMP$NEST.ID %in% NEST_wGen %>% sum) -7



NEST_wFam <- MORPHO %>% 
mutate( tmp = paste0(family1_morpho, family2_morpho, family1_visual, family2_visual),
		tmp2 = str_detect(tmp, 'f:')) %>% filter(tmp2 == TRUE) %>% pull(ID_nest) 

TOT.FAM <- (RES.COMP$NEST.ID %in% NEST_wFam %>% sum) -7

NEST_wOrd <- MORPHO %>% 
mutate( tmp = paste0(order1_morpho, order2_morpho, order1_visual, order2_visual),
		tmp2 = str_detect(tmp, 'o:'))  %>% filter(tmp2 == TRUE) %>% pull(ID_nest) 

TOT.ORD <- (RES.COMP$NEST.ID %in% NEST_wOrd %>% sum) -7



## plotting

comparison_plot <- function(data, x, y, label) {
  ggplot(data, aes({{x}}, {{y}})) +
    geom_col(position = position_dodge(), colour = 'black') +
	ylab('Samples with agreement [%]')+
	theme_bw() + theme(panel.grid.major=element_blank(),
		panel.grid.minor=element_blank(), axis.line=element_line(colour = 'black'))+
	geom_text(aes({{x}}, {{y}}, label = {{label}}), position = position_dodge(width = 1), vjust = -0.5)+
	theme(legend.position = "none") +
	geom_hline(yintercept = 50, col = "gray") 
		}

#across all samples
P.COMP1 <- RES.COMP2 %>% group_by(tax.level) %>%
		summarise(n.cor = sum(ID.cor)) %>%
		mutate(total = c(TOT.ORD, TOT.FAM, TOT.GEN, TOT.SP), pct = round((n.cor/total *100), digits =1)) %>%
			comparison_plot(., tax.level, pct, n.cor) +
			xlab('taxonomic level')+ scale_y_continuous(limits = c(0,110), breaks = c(0,25,50,100)) +
			annotate(geom = "text", x = c(1:4), y = 100, label = c(TOT.ORD,TOT.FAM,TOT.GEN,TOT.SP))


#experiment 1 only
P.COMP2 <- RES.COMP2 %>% filter(EXP %in% c('swab', 'feces')) %>% 
		filter(NEST.ID %in% c(NEST_wOrd, NEST_wFam, NEST_wGen, NEST_wSp)) %>%
		group_by(Sample.Type, NEST.Size, tax.level) %>%
		summarise(n.cor = sum(ID.cor),
			total = n()) %>%
		mutate(pct = round((n.cor/total *100), digits =1)) %>%
		comparison_plot(., NEST.Size, pct, n.cor) +
			facet_grid(Sample.Type ~ tax.level) +
			xlab('nest size') +
			geom_vline(xintercept = c(1.5, 2.5), lty = 'dashed', col = 'gray')+
			scale_x_discrete(labels=c("S", "M", "L")) +
			scale_y_continuous(limits = c(0,110), breaks = c(0,25,50,100)) 
			 


#experiment 2 only
P.COMP3 <- RES.COMP2 %>% filter(EXP %in% c('para', 'mixed')) %>% 
		filter(NEST.ID %in% c(NEST_wOrd, NEST_wFam, NEST_wGen, NEST_wSp)) %>%
		group_by(EXP, NEST.Size, tax.level) %>%
		summarise(n.cor = sum(ID.cor),
			total = n()) %>%
		mutate(pct = round((n.cor/total *100), digits =1)) %>%
		comparison_plot(., NEST.Size, pct, n.cor) +
			facet_grid(EXP ~ tax.level) +
		xlab('nest size')+
		geom_vline(xintercept = c(1.5, 2.5), lty = 'dashed', col = 'gray')+
		scale_x_discrete(labels=c("S", "M", "L"))+
			scale_y_continuous(limits = c(0,110), breaks = c(0,25,50,100)) 



plot_grid(P.COMP1, P.COMP2, P.COMP3, ncol =1, labels = c('A', 'B', 'C'))

bottom_row <- plot_grid(P.COMP2, P.COMP3, labels = c('B', 'C'))
plot_grid(P.COMP1, bottom_row, labels = c('A', ''), ncol = 1)


ggsave('HY.COMP.pdf', height = 10, width = 7)
ggsave('HY.COMP.jpg', height = 10, width = 7)









##SUPPLEMENTAL FIGURE S1
## DO WITH PREVIOS FUNCTION?

plot_suppl <- function(data, x, y, label, fill){
		ggplot(data, aes({{x}}, {{y}}, fill = {{fill}})) +
		geom_col(position = position_dodge(), colour = 'black', width = 0.75)+
		theme_bw() + theme(panel.grid.major=element_blank(),
		panel.grid.minor=element_blank(), axis.line=element_line(colour = 'black'))+
		xlab('No. of brood cells') + ylab('Samples with Hymenoptera \n at species level [%]')+
		geom_text(aes(x = {{x}}, y = {{y}}, label = {{label}}), position = position_dodge(width = 1), vjust = -0.5)+
		geom_hline(yintercept = 50)+
				geom_vline(xintercept = c(6, 12), lty = 'dashed', col = 'gray')+
		scale_y_continuous(breaks = c(0,25,50,75,100))+
		 expand_limits(x =0, y = 110)
		
		}





HY.SPEC.SIZE1 <- HY.tib %>% filter(Exp.Group %in% c('feces', 'swab')) %>% group_by(No.Cells, SampleType) %>%
	summarise(n.HySp = sum(Hy.det == 1),
		total = n()) %>%
	mutate(pct = round((n.HySp/total*100), digits = 1),
		pseudo.no = ifelse(No.Cells < 6, No.Cells, 
					ifelse(No.Cells < 11, No.Cells+1, No.Cells+2))	)
	



P.Hy.SIZE1 <- plot_suppl(HY.SPEC.SIZE1, pseudo.no, pct, n.HySp, SampleType) +
		facet_wrap(~SampleType) +
		theme(legend.position = "top")+
		scale_x_continuous(breaks = c(1:5, 7:11, 13:14, 17),
					labels = c(1:5, 6:10, 11, 12, 15))+
		theme(strip.background = element_blank(), 
		strip.text.x = element_blank()) +
		scale_fill_discrete(name = 'Sample type', labels = c('feces', 'swabs'))+
		annotate(geom = "text", x = 3, y = 110, label = 'N = 8')+
		annotate(geom = "text", x = 9, y = 110, label = 'N = 8')+
		annotate(geom = "text", x = 15, y = 110, label = 'N = 10 / 8 / 2')
	



HY.SPEC.SIZE2 <- HY.tib %>% filter(Exp.Group == 'para') %>% group_by(No.Cells, SampleType) %>%
	summarise(n.HySp = sum(Hy.det == 1),
		total = n()) %>%
	mutate(pct = round((n.HySp/total*100), digits = 1),
		pseudo.no = ifelse(No.Cells < 6, No.Cells, 
					ifelse(No.Cells < 11, No.Cells+1, No.Cells+2))	)
	


P.Hy.SIZE2 <- plot_suppl(HY.SPEC.SIZE2, pseudo.no, pct, n.HySp, SampleType) +
		theme(legend.position = "none")+
		scale_x_continuous(breaks = c(1:5, 7:11, 13:14),
					labels = c(1:5, 6:10, 11, 12))+
		theme(strip.background = element_blank(), 
		strip.text.x = element_blank()) +
		scale_fill_manual(name = 'Sample type', labels = c('feces', 'swabs'), values = 'grey')+
		annotate(geom = "text", x = 3, y = 110, label = 'N = 4 / 6 / 4 / 20 / 8')+
		annotate(geom = "text", x = 9, y = 110, label = 'N = 4 / 6 / 2 / 4 / 8')+
		annotate(geom = "text", x = 14, y = 110, label = 'N = 4 / 4')

	


HY.SPEC.SIZE3 <- HY.tib %>% filter(Exp.Group == 'mixed') %>% group_by(No.Cells, SampleType) %>%
	summarise(n.HySp = sum(Hy.det == 1),
		total = n()) %>%
	mutate(pct = round((n.HySp/total*100), digits = 1),
		pseudo.no = ifelse(No.Cells < 6, No.Cells, 
					ifelse(No.Cells < 11, No.Cells+1, No.Cells+2))	)
	


P.Hy.SIZE3 <- plot_suppl(HY.SPEC.SIZE3, pseudo.no, pct, n.HySp, SampleType) +
		theme(legend.position = "none")+
		scale_x_continuous(breaks = c(3:5, 7:8, 13, 17),
 					labels = c(3:5, 6:7, 11, 15))+
		theme(strip.background = element_blank(), 
		strip.text.x = element_blank()) +
		scale_fill_manual(name = 'Sample type', labels = c('feces', 'swabs'), values = 'grey')+
		annotate(geom = "text", x = 3, y = 110, label = 'N = 14 / 4 / 8')+
		annotate(geom = "text", x = 9, y = 110, label = 'N = 4 / 2 ')+
		annotate(geom = "text", x = 14, y = 110, label = 'N = 2 / 2')



plots <- align_plots(P.Hy.SIZE1, P.Hy.SIZE2, align = 'v', axis = '1')
bottom_row <- plot_grid(plots[[2]], P.Hy.SIZE3, labels = c('B', 'C'))

plot_grid(plots[[1]], bottom_row, labels = c('A', ''), ncol = 1)

ggsave('HY.DET.size.pdf', height = 7.5, width = 10)
ggsave('HY.DET.size.jpg', height = 7.5, width = 10)


