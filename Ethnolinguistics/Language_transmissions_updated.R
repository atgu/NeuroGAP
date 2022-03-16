## looking into the relationship between genetics and language in NeuroGAP


#load in needed packages 
## Allan - note you will also need to load the dplyr library for file 

#install.packages("cowplot")
#install.packages("dplyr")
#install.packages("plyr")
#install.packages("ggplot2")
#install.packages("naniar")
#you will only need to install the package once per R environment. You'll still need to load the library each time though

library(ggplot2)
library(cowplot)
library(dplyr) 
library(plyr)
library(naniar)


setwd('/Users/elizabeth/OneDrive - Baylor College of Medicine/Projects/MGH-Broad/NeuroGap/GeneticsLanguage')
setwd('/Users/HP/Desktop/LANGUAGE_GINGER/')

# Load in language phenotypes
## Here you will also need to read in the sample based file (neurolang)
langphenos <- read.csv('LanguageCodes_withLangFamilies-1.csv', header = T)
neuroLang <- read.csv('LanguageEthnicity_phenos_NeuroGAP_v2.csv', header = T)

##Correction
neuroLang2 <- read.csv('NGP_Freeze3_EthnLangUpdate.csv', header = T)  

# Joining the sample level data with Yakov's language info
## Here you are taking langphenos (data from Yakov) and joining it with neuroLang (sample level data) 
## By doing a left join you are keeping everything in the neurolang file and adding on the langpheno data 
consent_lang <- neuroLang %>% left_join(langphenos, by = c("consent_lang" = "Code"))

#write out a table of just the shareable info for the manuscript
ms_langphenos <- neuroLang[,-c(1)]



##Correction
consent_lang2 <- neuroLang2 %>% right_join(langphenos, by = c("consent_lang" = "Code"))


#plot up the newly assigned language overall. This assumes there's a field called 'Langauge' that we will be wanting to color bars by
p <- ggplot(consent_lang, aes(x=Langauge, fill=Langauge)) + ggtitle("Assigned language") + xlab("Assigned language") + ylab("Count") + geom_histogram(stat="count", position = "stack") +  theme_bw()
p #take a look at the plot in your R studio screen
# The command below will save the plot in your current working directory, which you set at the beginning
# you can play around with how big this needs to be to look good when exported by changing the height and width
save_plot('NeuroGAP_Assigned_lang_overall.pdf', p, base_height=7, base_width=14)

#plot up the newly assigned language by study country. Assuming field name is 'study_country'
p1 <- ggplot(consent_lang, aes(x=Langauge, fill=study_country)) + ggtitle("Assigned language") + xlab("Assigned Language") + ylab("Count") + geom_histogram(stat="count", position = "stack") +  theme_bw() 
p1 #take a look at the plot in your R studio screen
save_plot('NeuroGAP_Assigned_lang_country.pdf', p1, base_height=7, base_width=14) 

#plot up the the four countries filled by the newly assigned languages 
p2 <- ggplot(consent_lang, aes(x=study_country, fill=Langauge)) + ggtitle("Study Countries Assigned Language") + xlab("Study Country") + ylab("Count") + geom_histogram(stat="count", position = "stack") +  theme_bw() 
p2 #take a look at the plot in your R studio screen
save_plot('NeuroGAP_Study_Country_languages.pdf', p2, base_height=7, base_width=14) 

p3 <- ggplot(consent_lang, aes(x=study_country, fill=lang_family_manual)) + ggtitle("Study Countries Manually Assigned Language") + xlab("Study Country") + ylab("Count") + geom_histogram(stat="count", position = "stack") +  theme_bw() 
p3


## Plots excluding English

Px <- ggplot (subset(consent_lang, Language %in% c("English")), aes(x=Language, fill=Language)) + ggtitle("Assigned language") + xlab("Assigned language") + ylab("Count") + geom_histogram(stat="count", position = "stack") +  theme_bw()


No_English <- consent_lang[!grepl("English", consent_lang$Langauge),]

## Plot: Newly assigned language overall, excluding English

P4 <- ggplot(No_English, aes(x=Langauge, fill=Langauge)) + ggtitle("Assigned language excluding English") + xlab("Assigned language") + ylab("Count") + geom_histogram(stat="count", position = "stack") +  theme_bw()
P4 #take a look at the plot in your R studio screen
# The command below will save the plot in your current working directory, which you set at the beginning
save_plot('Assigned_lang_overall_Excluding English.pdf', p4, base_height=7, base_width=14)

## Plot: Newly assigned language by study country. Assuming field name is 'study_country'

p5 <- ggplot(No_English, aes(x=Langauge, fill=study_country)) + ggtitle("Assigned language by Country excluding English") + xlab("Assigned Language") + ylab("Count") + geom_histogram(stat="count", position = "stack") +  theme_bw() 
p5 #take a look at the plot in your R studio screen
save_plot('NeuroGAP_Assigned_lang_country.pdf', p5, base_height=7, base_width=14) 

## Plot: The four countries filled by the newly assigned languages, excluding English 

p6 <- ggplot(No_English, aes(x=study_country, fill=Langauge)) + ggtitle("Study Countries Assigned Language excluding English") + xlab("Study Country") + ylab("Count") + geom_histogram(stat="count", position = "stack") +  theme_bw() 
p6 #take a look at the plot in your R studio screen
save_plot('NeuroGAP_Study_Country_languages.pdf', p6, base_height=7, base_width=14)

## Plot: Study countries manually assigned language excluding English

p7 <- ggplot(No_English, aes(x=study_country, fill=lang_family_manual)) + ggtitle("Study Countries Manually Assigned Language excluding English") + xlab("Study Country") + ylab("Count") + geom_histogram(stat="count", position = "stack") +  theme_bw() 


############################################# Language transmission #################################################

# Paternal Language transmission

length(which(na.exclude(consent_lang$lang_self_1) ==consent_lang$lang_pat_1))
# 2562


#how many entries in lang_pat_1 are not NA

length(na.exclude(consent_lang$lang_pat_1))
#3675

#how many IDs are not NA for both lang_self_1 and lang_pat_1

langpat = consent_lang[c("lang_self_1", "lang_pat_1")]

dim(na.exclude(langpat))
pat_trans_rate <- length(which(na.exclude(consent_lang$lang_self_1) ==consent_lang$lang_pat_1))/length(na.exclude(consent_lang$lang_pat_1))
#0.6971429

## now re-do but with downsampling
pat_trans_rate <- length(which(na.exclude(consent_lang$lang_self_1) ==consent_lang$lang_pat_1))/length(na.exclude(consent_lang$lang_pat_1))



########################################## Maternal language transmission ########################################

## Maternal transmission including English
# Maternal language transmission

length(which(na.exclude(consent_lang$lang_self_1) ==consent_lang$lang_mat_1))
#2324

#how many entries in lang_mat_1 are not NA

length(na.exclude(consent_lang$lang_mat_1))
#3781

#how many IDs are not NA for both lang_self_1 and lang_mat_1

langmat = consent_lang[c("lang_self_1", "lang_mat_1")]

dim(na.exclude(langmat))

mat_trans_rate <- length(which(na.exclude(consent_lang$lang_self_1) ==consent_lang$lang_mat_1))/length(na.exclude(consent_lang$lang_mat_1))
#0.6146522


########################################## Paternal grand mothers language transmission ########################################
## Paternal grand mothers language transmission

length(which(na.exclude(consent_lang$lang_self_1) ==consent_lang$lang_pgm_1))


#how many entries in lang_pm_1 are not NA

length(na.exclude(consent_lang$lang_pgm_1))


#how many IDs are not NA for both lang_self_1 and lang_pgm_1

langpgm = consent_lang[c("lang_self_1", "lang_pgm_1")]

dim(na.exclude(langpgm))

########################################## Paternal grand fathers language transmission ########################################

length(which(na.exclude(consent_lang$lang_self_1) ==consent_lang$lang_pgf_1))


#how many entries in lang_pgf_1 are not NA

length(na.exclude(consent_lang$lang_pgf_1))


#how many IDs are not NA for both lang_self_1 and lang_pgf_1

langpgf = consent_lang[c("lang_self_1", "lang_pgf_1")]

dim(na.exclude(langpgf))


########################################## Maternal grand mothers language transmission ########################################

length(which(na.exclude(consent_lang$lang_self_1) ==consent_lang$lang_mgm_1))


#how many entries in lang_mgm_1 are not NA

length(na.exclude(consent_lang$lang_mgm_1))


#how many IDs are not NA for both lang_self_1 and lang_mgm_1

langmgm = consent_lang[c("lang_self_1", "lang_mgm_1")]

dim(na.exclude(langmgm))

########################################## Maternal grand mothers language transmission ########################################

length(which(na.exclude(consent_lang$lang_self_1) ==consent_lang$lang_mgf_1))


#how many entries in lang_mgf_1 are not NA

length(na.exclude(consent_lang$lang_mgf_1))


#how many IDs are not NA for both lang_self_1 and lang_mgf_1

langmgf = consent_lang[c("lang_self_1", "lang_mgf_1")]

dim(na.exclude(langmgf))




######################## Partitioning language transmissions by Matrilineal vs Patrillineal ######################

#### UPDATE for reviewss - correction for sample size
# N=105 and 674 for matrilineal and patrilineal transmission respectively
#keep matrilineal as is. Downsample the patrilineal to be the same length and re-do transmissions to make sure things are still consistent



# load cultural information vector

cultural <- read.csv('inheritance_EA076.csv', header= T, sep= ",")
#cultural <- read.csv('/Users/elizabeth/OneDrive - Baylor College of Medicine/Projects/MGH-Broad/NeuroGap/GeneticsLanguage/NeuroGAP_GeneticsLanguage/Results/Language_transmission_plots/inheritance_EA076.csv', header= T, sep= ",")


## Merging the  files: By doing a left join we are keeping everything in the inheritance.csv file and adding on the consent)lang data 

merged_maternal_paternal <- merge(consent_lang, cultural, by = "IID")


#################### partitioning to matrilineal ##################

matrilineal <- merged_maternal_paternal [which (merged_maternal_paternal$inheritance_EA076 =="Matrilineal"),]

################### partitioning to patrilineal ##################

patrilineal <- merged_maternal_paternal [which (merged_maternal_paternal$inheritance_EA076 =="Patrilineal"),]

#randomly grab the same number of rows as is matrilineal - 105. We can use the sample_n command for this. 
patri_subsample <- sample_n(patrilineal, size = 105, replace= F, prob=NULL)
#Then we can re-do the patri analyses with patri_subsample instead of patrilineal as the dataset to have a matched sample!




# Follow previous steps for language transmissions.

############################################# Language transmission: Patrilineal #################################################

# Paternal Language transmission

length(which(na.exclude(patrilineal$lang_self_1) ==patrilineal$lang_pat_1))


#how many entries in lang_pat_1 are not NA

length(na.exclude(patrilineal$lang_pat_1))


#how many IDs are not NA for both lang_self_1 and lang_pat_1

langpat = patrilineal[c("lang_self_1", "lang_pat_1")]

dim(na.exclude(langpat))


########################################## Maternal language transmission ########################################

## Maternal transmission including English
# Maternal language transmission

length(which(na.exclude(patrilineal$lang_self_1) ==patrilineal$lang_mat_1))


#how many entries in lang_mat_1 are not NA

length(na.exclude(patrilineal$lang_mat_1))


#how many IDs are not NA for both lang_self_1 and lang_mat_1

langmat = patrilineal[c("lang_self_1", "lang_mat_1")]

dim(na.exclude(langmat))


########################################## Paternal grand mothers language transmission ########################################
## Paternal grand mothers language transmission

length(which(na.exclude(patrilineal$lang_self_1) ==patrilineal$lang_pgm_1))


#how many entries in lang_pm_1 are not NA

length(na.exclude(patrilineal$lang_pgm_1))


#how many IDs are not NA for both lang_self_1 and lang_pgm_1

langpgm = patrilineal[c("lang_self_1", "lang_pgm_1")]

dim(na.exclude(langpgm))

########################################## Paternal grand fathers language transmission ########################################

length(which(na.exclude(patrilineal$lang_self_1) ==patrilineal$lang_pgf_1))


#how many entries in lang_pgf_1 are not NA

length(na.exclude(patrilineal$lang_pgf_1))


#how many IDs are not NA for both lang_self_1 and lang_pgf_1

langpgf = patrilineal[c("lang_self_1", "lang_pgf_1")]

dim(na.exclude(langpgf))


########################################## Maternal grand mothers language transmission ########################################

length(which(na.exclude(patrilineal$lang_self_1) ==patrilineal$lang_mgm_1))


#how many entries in lang_mgm_1 are not NA

length(na.exclude(patrilineal$lang_mgm_1))


#how many IDs are not NA for both lang_self_1 and lang_mgm_1

langmgm = patrilineal[c("lang_self_1", "lang_mgm_1")]

dim(na.exclude(langmgm))

########################################## Maternal grand mothers language transmission ########################################

length(which(na.exclude(patrilineal$lang_self_1) ==patrilineal$lang_mgf_1))


#how many entries in lang_mgf_1 are not NA

length(na.exclude(patrilineal$lang_mgf_1))


#how many IDs are not NA for both lang_self_1 and lang_mgf_1

langmgf = patrilineal[c("lang_self_1", "lang_mgf_1")]

dim(na.exclude(langmgf))

############################################# Language transmission: Matrilineal #################################################

# Paternal Language transmission

length(which(na.exclude(matrilineal$lang_self_1) ==matrilineal$lang_pat_1))


#how many entries in lang_pat_1 are not NA

length(na.exclude(matrilineal$lang_pat_1))


#how many IDs are not NA for both lang_self_1 and lang_pat_1

langpat = matrilineal[c("lang_self_1", "lang_pat_1")]

dim(na.exclude(langpat))


########################################## Maternal language transmission ########################################

## Maternal transmission including English
# Maternal language transmission

length(which(na.exclude(matrilineal$lang_self_1) ==matrilineal$lang_mat_1))


#how many entries in lang_mat_1 are not NA

length(na.exclude(matrilineal$lang_mat_1))


#how many IDs are not NA for both lang_self_1 and lang_mat_1

langmat = matrilineal[c("lang_self_1", "lang_mat_1")]

dim(na.exclude(langmat))


########################################## Paternal grand mothers language transmission ########################################
## Paternal grand mothers language transmission

length(which(na.exclude(matrilineal$lang_self_1) ==matrilineal$lang_pgm_1))


#how many entries in lang_pm_1 are not NA

length(na.exclude(matrilineal$lang_pgm_1))


#how many IDs are not NA for both lang_self_1 and lang_pgm_1

langpgm = matrilineal[c("lang_self_1", "lang_pgm_1")]

dim(na.exclude(langpgm))

########################################## Paternal grand fathers language transmission ########################################

length(which(na.exclude(matrilineal$lang_self_1) ==matrilineal$lang_pgf_1))


#how many entries in lang_pgf_1 are not NA

length(na.exclude(matrilineal$lang_pgf_1))


#how many IDs are not NA for both lang_self_1 and lang_pgf_1

langpgf = matrilineal[c("lang_self_1", "lang_pgf_1")]

dim(na.exclude(langpgf))


########################################## Maternal grand mothers language transmission ########################################

length(which(na.exclude(matrilineal$lang_self_1) ==matrilineal$lang_mgm_1))


#how many entries in lang_mgm_1 are not NA

length(na.exclude(matrilineal$lang_mgm_1))


#how many IDs are not NA for both lang_self_1 and lang_mgm_1

langmgm = matrilineal[c("lang_self_1", "lang_mgm_1")]

dim(na.exclude(langmgm))

########################################## Maternal grand mothers language transmission ########################################

length(which(na.exclude(patrilineal$lang_self_1) ==matrilineal$lang_mgf_1))


#how many entries in lang_mgf_1 are not NA

length(na.exclude(matrilineal$lang_mgf_1))


#how many IDs are not NA for both lang_self_1 and lang_mgf_1

langmgf = matrilineal[c("lang_self_1", "lang_mgf_1")]

dim(na.exclude(langmgf))


###########Corrected language transmissions

############################################# Language transmission_corrected #################################################

# Corrected Paternal Language transmission
consent_lang3 <- consent_lang2 %>% dplyr::na_if(-888)

length(which(na.exclude(consent_lang3$lang_self_1) ==consent_lang3$lang_pat_1))


#how many entries in lang_pat_1 are not NA

length(na.exclude(consent_lang3$lang_pat_1))


#how many IDs are not NA for both lang_self_1 and lang_pat_1

langpat3 = consent_lang3[c("lang_self_1", "lang_pat_1")]

dim(na.exclude(langpat))


########################################## Maternal language transmission ########################################

## Maternal transmission including English
# Maternal language transmission

length(which(na.exclude(consent_lang$lang_self_1) ==consent_lang$lang_mat_1))


#how many entries in lang_mat_1 are not NA

length(na.exclude(consent_lang$lang_mat_1))


#how many IDs are not NA for both lang_self_1 and lang_mat_1

langmat = consent_lang[c("lang_self_1", "lang_mat_1")]

dim(na.exclude(langmat))


########################################## Paternal grand mothers language transmission ########################################
## Paternal grand mothers language transmission

length(which(na.exclude(consent_lang$lang_self_1) ==consent_lang$lang_pgm_1))


#how many entries in lang_pm_1 are not NA

length(na.exclude(consent_lang$lang_pgm_1))


#how many IDs are not NA for both lang_self_1 and lang_pgm_1

langpgm = consent_lang[c("lang_self_1", "lang_pgm_1")]

dim(na.exclude(langpgm))

########################################## Paternal grand fathers language transmission ########################################

length(which(na.exclude(consent_lang$lang_self_1) ==consent_lang$lang_pgf_1))


#how many entries in lang_pgf_1 are not NA

length(na.exclude(consent_lang$lang_pgf_1))


#how many IDs are not NA for both lang_self_1 and lang_pgf_1

langpgf = consent_lang[c("lang_self_1", "lang_pgf_1")]

dim(na.exclude(langpgf))


########################################## Maternal grand mothers language transmission ########################################

length(which(na.exclude(consent_lang$lang_self_1) ==consent_lang$lang_mgm_1))


#how many entries in lang_mgm_1 are not NA

length(na.exclude(consent_lang$lang_mgm_1))


#how many IDs are not NA for both lang_self_1 and lang_mgm_1

langmgm = consent_lang[c("lang_self_1", "lang_mgm_1")]

dim(na.exclude(langmgm))

########################################## Maternal grand mothers language transmission ########################################

length(which(na.exclude(consent_lang$lang_self_1) ==consent_lang$lang_mgf_1))


#how many entries in lang_mgf_1 are not NA

length(na.exclude(consent_lang$lang_mgf_1))


#how many IDs are not NA for both lang_self_1 and lang_mgf_1

langmgf = consent_lang[c("lang_self_1", "lang_mgf_1")]

dim(na.exclude(langmgf))

## Here are some tips for changing and manipulating what is being plotted. You may know some of this already but I thought it would be helful to include.
# If you would like to edit what is being plotted, you can change the 'x' and 'fill' arguments in aes() in the ggplot() function. 
# As is suggested by the name, 'x' will be what is shown on the x axis, and fill is what the histogram bars will be coloured in with. 
# You can change these to any of the columns in the consent_lang variable. To see what is in consent_lang you can click on 'consent_lang'
# in the environment tab under Data in the upper right corner and that will open a tab in the main window so you can see the data.  
# In regards to joining the two files, I think the current command will allow you to plot the data in the format you were suggesting 
# I included a 3rd plot (p2) command which I think is doing what you had suggested, but if you meant something different, let me know. 

