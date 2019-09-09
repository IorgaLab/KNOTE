

#Be sure to load all given libraries and functions in the first categories. Data loading is used to import the dataset in the file.
#Make sure to set the working directory where the dataset file and the function script are put in your system, an example is available below

#setwd("C:/Users/[Username]/PATH/TO/FILE")

###### Libraries ######
library(caret)
library(randomForest)
library(ggplot2)
library(ROCR)

###### Functions ######
# Removing the individuals for which MICs of the selected antibiotic were determined
rem_ND <- function(data_toclean, resist_column){
    if (length(which(data_toclean[,resist_column] == "ND")) != 0){
        data_cleaned = data_toclean[-which(data_toclean[,resist_column] == "ND"),]
    }
    else{
        data_cleaned = data_toclean
    }
    return(data_cleaned)
}
###### Data loading ######
#The data must exist as a .csv file to be correctly read by the program
data = read.csv("3b_390genomes.csv", header = T, sep = ";", stringsAsFactors = F)
data = data[-which(data[,7976:7987] == "ND"),]

#Names for every resistance column in the original dataset
colnames(data)[7984] = "Meropenem_SR"
colnames(data)[7985] = "Levofloxacin_SR"
colnames(data)[7986] = "Amikacin_SR"
colnames(data)[7987] = "Colistin_SR"
colnames(data)[7980] = "Meropenem_let"
colnames(data)[7981] = "Levofloxacin_let"
colnames(data)[7982] = "Amikacin_let"
colnames(data)[7983] = "Colistin_let"
colnames(data)[7976] = "Meropenem_log2"
colnames(data)[7977] = "Levofloxacin_log2"
colnames(data)[7978] = "Amikacin_log2"
colnames(data)[7979] = "Colistin_log2"

colnames(data)[1] = "GenomeCodes"

listColNames = colnames(data[,which(regmatches(colnames(data), regexpr(".$", colnames(data))) == ".")])
listNumColNames = which(regmatches(colnames(data), regexpr(".$", colnames(data))) == ".")
for (name in listColNames){
    name = paste(substr(name, 1, nchar(name)-1), "DEL", sep = "")
    for (number in listNumColNames){
        colnames(data)[number] = name
    }
}



###### Enrich ######

PrepEnrich <- function(data, column, threshold){
##Function aiming at creating sub-datasets from an original dataset of resistance motifs and continuous resistance.
##It will create intervals containing the closest values 2 by 2 and removing all common motifs between all individuals from both values.
##It creates a list of datasets for each interval that will be used to oversample the resistance values that are unrepresented in the
##original dataset.
    
##Data : Dataset that will be separated
##Column : The resistance column
##Threshold : Percentage (between 0 and 1) under which individuals for a given resistance value will be considered undersampled from the
##most represented value in the dataset
    #Creating a graph to see which columns will or can be selected
    data = rem_ND(data, column)
    data = cbind(data[,1:7975], data[,column])
    codes = data[,1]
    colnames(data)[7976] = "Antibiotic"
    hist_antibiotic = ggplot(data=data, aes(reorder(data$Antibiotic))) + stat_count(width = 1) + xlab("MIC") + ylab("Number of individuals")
    percentage1 = 0.10*max(table(data[,7976]))
    hist_antibiotic = hist_antibiotic + geom_hline(yintercept=percentage1, linetype="dashed", color = "orange") + geom_text(aes(1, percentage1,label = "10%"), vjust = -1, color = "orange")
    percentage2 = 0.25*max(table(data[,7976]))
    hist_antibiotic = hist_antibiotic + geom_hline(yintercept=percentage2, linetype="dashed", color = "orange") + geom_text(aes(1, percentage2,label = "25%"), vjust = -1, color = "orange")
    percentage3 = 0.50*max(table(data[,7976]))
    hist_antibiotic = hist_antibiotic + geom_hline(yintercept=percentage3, linetype="dashed", color = "orange") + geom_text(aes(1, percentage3,label = "50%"), vjust = -1, color = "orange")
    percentage4 = 0.75*max(table(data[,7976]))
    hist_antibiotic = hist_antibiotic + geom_hline(yintercept=percentage4, linetype="dashed", color = "orange") + geom_text(aes(1, percentage4,label = "75%"), vjust = -1, color = "orange")
    percthresh = threshold*max(table(data[,7976]))
    labpercentage = paste(as.character(threshold*100), "%", sep = "")
    hist_antibiotic = hist_antibiotic + geom_hline(yintercept=percthresh, linetype="dashed", color = "red") + geom_text(aes(1, percthresh,label = labpercentage), vjust = -1, color = "red")
    plot(hist_antibiotic)
        
    cmi = sort(as.numeric(names(table(data[,7976]))))
    ##Creating the list of intervals
    intervals = list()
    index = 1
    for(i in seq(1:(length(cmi)-1))){
        infbound = sum(data[,7976] == cmi[i])
        supbound = sum(data[,7976] == cmi[i+1])
        if ((infbound < threshold*max(table(data[,7976]))) || (supbound < threshold*max(table(data[,7976])))){
            temp = c(cmi[i], cmi[i+1])
            intervals[[index]] = temp
            index = index + 1
            cat("Added : [", cmi[i], ",", cmi[i+1], "] to the list of intervals \n")
        }
        
    }
    ##Counting and creating the lists of individuals for each interval
    AllSubsets = list()
    for (i in seq(1:(length(intervals)))){
        values_inf = which(data[,7976] == intervals[[i]][1])
        values_sup = which(data[,7976] == intervals[[i]][2])
        AllSubsets[[i]] = rbind(data[values_inf,], data[values_sup,])
        cat("There are", nrow(AllSubsets[[i]]), "individuals in the [", intervals[[i]], "] interval \n")
    }
    table(data[,7976])
    
    ##Removing the common descriptors between every individual in each sub-dataset
    for (i in seq(1:(length(AllSubsets)))){
        toremove = numeric()
        for (j in seq(2:7976)){
            if (length(unique(AllSubsets[[i]][,j])) < 2){
                toremove = c(toremove, j)
            }
        }
        cat("Removed", length(toremove), "descriptors from the [", intervals[[i]], "] interval \n")
        AllSubsets[[i]] = AllSubsets[[i]][,-toremove]
    }
    return(AllSubsets) ##Return the list of datasets for further use
}    

Enrich <- function(data_target, data_codes, AllSubsets){
    ##This function aims at enriching every subset created earlier so that the final dataset will be larger and more useful in predictions
    ##data_target : The final dataset that will contain every new individual
    ##data_codes : The original genome codes
    ##AllSubsets : List of sub-datasets created by the PrepEnrich function earlier
    
    for (value in seq(1:length(AllSubsets))){ 
        ##For every sub-dataset, we will isolate each one of them and treat them separately
        #Isolation of a sub-dataset, named currDF (for "current DataFrame")
        currDF = AllSubsets[[value]][,-1] 
        currDFcodes = AllSubsets[[value]][,1] ##Genome codes for the CurrDF
        ncolcurrDF = ncol(currDF)
        currDF[,ncolcurrDF] = droplevels(currDF[,ncolcurrDF])
        
        #Calculation of the number of individuals in each value of the subset, to determine which will be enriched
        uniqcurrDF = sort(as.numeric(as.character(unique(currDF[,ncol(currDF)]))))
        len_leftinter = length(which(data_target[,7976] == uniqcurrDF[1]))
        print(len_leftinter)
        len_rightinter = length(which(data_target[,7976] == uniqcurrDF[2]))
        print(len_rightinter)
        
        #We realise a random forest prediction to determine the impactful motifs in the difference between the two values
        #so that they are not changed in the creation of new individuals
        RFmino = randomForest(Antibiotic ~ ., data = currDF, ntree = 20, mtry = 10, importance = T)   
        currDF.importance = abs(RFmino$importance[,3])
        currDF.importance = sort(currDF.importance, decreasing = T)
        #The "toadd" variable contains how much individuals must be created for :
        #Either having the same number of the most represented resistance value in the original dataset
        #Or the highest number of unique individuals that can be created from the unimportant motifs determined by the previous 
        #random forest
        toadd = min(max(table(data_target[,ncol(data_target)])) - min(len_rightinter, len_leftinter), length(which(currDF.importance == 0)))
        if (toadd < max(table(data_target[,ncol(data_target)])) - min(len_rightinter, len_leftinter)){
            cat("Length of descriptors was inferior to number of individuals needed, making as much individuals as possible from importance test \n")
        }
        tailimpt = tail(currDF.importance, toadd)
        
        #Selecting which interval will be enriched
        if (len_leftinter < len_rightinter){
            #Whether the left or right is enriched, the method is exactly the same :
            # Step 1) A random and original individual in the selected interval is chosen and it's information is kept in variables
            # Step 2) A random unimportant motif (in the tailimpt variable) is selected and removed 
            #         from the list (so that every new individual is unique too)
            # Step 3) That motif is put to 0 or 1, whether it was 1 or 0 in the original individual
            # Step 4) The presence or absence of the Wild-Type variation of the motif is checked if a mutated motif was selected
            #         (to prevent Wild-type and mutated motif to co-exist)
            # Step 5) The new individual is added to data_target, the final dataset, with exactly the same motifs as its original, except
            #        the one selected
            print("Enriching left interval")
            for (i in seq(1:toadd)){ #Step 1
                if (length(which(currDF[,ncolcurrDF] == uniqcurrDF[2])) > 1){
                    codeorig = as.character(currDFcodes[sample(which(currDF[,ncolcurrDF] == uniqcurrDF[1]), 1)])
                }
                else{
                    codeorig = currDFcodes[which(currDF[,ncolcurrDF] == uniqcurrDF[1])]
                }
                cat("Using ", as.character(codeorig), " as original individual \n")
                newind_code = which(data_codes == codeorig)
                newind = data_target[newind_code,]
                tomutate = names(tailimpt[i]) #Step 2
                print(tomutate)
                tomutateorig = which(colnames(data_target) == tomutate)
                if (newind[,tomutateorig] == 0){ #Step 3
                    newind[,tomutateorig] = 1
                    if(endsWith(tomutate, "ABSENT") == FALSE){ #Step 4
                        expreg = regexpr(".*\\.", tomutate, perl=TRUE) #Find the resistance motif used for the new individual
                        motif = regmatches(tomutate, expreg)
                        expregWT = regexpr(".*\\.ABSENT", colnames(newind[,grep(motif,colnames(newind))])) #Find the wild type motif of the selected motif
                        motifWT = regmatches(colnames(newind[,grep(motif,colnames(newind))]), expregWT) 
                        if(length(motifWT) != 0){
                            if (tomutate != motifWT){
                                motifWTorig = which(colnames(newind) == motifWT)
                                if (newind[,motifWTorig] == 1){
                                    newind[,motifWTorig] = 0
                                    print("Removed wild-type motif because of mutation")
                                }
                                else{
                                    print("No wildtype was found in original individual")
                                }
                            }
                        }
                    }
                    else{
                        print("Added Wild-type motif to sequence")
                    }
                }
                else{ #Step 3
                    newind[,tomutateorig] = 0
                    if(endsWith(tomutate, "ABSENT") == FALSE){ #Step 4
                        expreg = regexpr(".*\\.", tomutate, perl=TRUE) #Find the resistance motif used for the new individual
                        motif = regmatches(tomutate, expreg)
                        expregWT = regexpr(".*\\.ABSENT", colnames(newind[,grep(motif,colnames(newind))])) #Find the wild type motif of the selected motif
                        motifWT = regmatches(colnames(newind[,grep(motif,colnames(newind))]), expregWT) 
                        if(length(motifWT) != 0){
                            if (tomutate != motifWT){
                                motifWTorig = which(colnames(newind) == motifWT)
                                if (newind[,motifWTorig] == 0){ 
                                    newind[,motifWTorig] = 1
                                    print("Removed wild-type motif because of mutation")
                                }
                                else{
                                    print("No wildtype was found in original individual")
                                }
                            }
                        }
                    }
                    else{
                        print("Removed Wild-type motif from sequence")
                    }
                }
                newind[,1] = paste(codeorig, "-added-", i) #Step 5
                data_target = rbind(data_target, newind)
            }
        }
        else{ 
            print("Enriching right interval")
            for (i in seq(1:toadd)){ #Step 1
                if (length(which(currDF[,ncolcurrDF] == uniqcurrDF[2])) > 1){
                    codeorig = as.character(currDFcodes[sample(which(currDF[,ncolcurrDF] == uniqcurrDF[2]), 1)])
                }
                else{
                    codeorig = currDFcodes[which(currDF[,ncolcurrDF] == uniqcurrDF[2])]
                }
                cat("Using", as.character(codeorig), "as original individual \n")
                newind_code = which(data_codes == codeorig)
                newind = data_target[newind_code,]
                tomutate = names(tailimpt[i]) #Step 2
                print(tomutate)
                tomutateorig = which(colnames(data_target) == tomutate)
                if (newind[,tomutateorig] == 0){ #Step 3
                    newind[,tomutateorig] = 1
                    if(endsWith(tomutate, "ABSENT") == FALSE){ #Step 4
                        expreg = regexpr(".*\\.", tomutate, perl=TRUE) #Find the resistance motif used for the new individual
                        motif = regmatches(tomutate, expreg)
                        expregWT = regexpr(".*\\.ABSENT", colnames(newind[,grep(motif,colnames(newind))])) #Find the wild type motif of the selected motif
                        motifWT = regmatches(colnames(newind[,grep(motif,colnames(newind))]), expregWT) 
                        if(length(motifWT) != 0){
                            if (tomutate != motifWT){
                                motifWTorig = which(colnames(newind) == motifWT)
                                if (newind[,motifWTorig] == 1){
                                    newind[,motifWTorig] = 0
                                    print("Removed wild-type motif because of mutation")
                                }
                                else{
                                    print("No wildtype was found in original individual")
                                }
                            }
                        }
                    }
                    else{
                        print("Added Wild-type motif to sequence")
                    }
                }
                else{ #Step 3
                    newind[,tomutateorig] = 0
                    if(endsWith(tomutate, "ABSENT") == FALSE){ #Step 4
                        expreg = regexpr(".*\\.", tomutate, perl=TRUE) #Find the resistance motif used for the new individual
                        motif = regmatches(tomutate, expreg)
                        expregWT = regexpr(".*\\.ABSENT", colnames(newind[,grep(motif,colnames(newind))])) #Find the wild type motif of the selected motif
                        motifWT = regmatches(colnames(newind[,grep(motif,colnames(newind))]), expregWT) 
                        if(length(motifWT) != 0){
                            if (tomutate != motifWT){
                                motifWTorig = which(colnames(newind) == motifWT)
                                if (newind[,motifWTorig] == 0){
                                    newind[,motifWTorig] = 1
                                    print("Removed wild-type motif because of mutation")
                                }
                                else{
                                    print("No wildtype was found in original individual")
                                }
                            }
                        }
                    }
                    else{
                        print("Removed Wild-type motif from sequence")
                    }
                }
                newind[,1] = paste(codeorig, "-added-", i) #Step 5
                data_target = rbind(data_target, newind)
            }
        }
    }
    colnames(data_target)[7976] = "Antibiotic" # A new histogram is created to compare with the first created in PrepEnrich
    hist_antibiotic = ggplot(data=data_target, aes(reorder(data_target$Antibiotic))) + stat_count(width = 1) + xlab("MIC") + ylab("Number of individuals")
    plot(hist_antibiotic)
    return(data_target) # The final dataset is returned
}


AllSubsets = PrepEnrich(data, 7979, 0.1) ##Calling the first function

data_target = cbind(data[,1:7975], data[,7979]) ##Creating the soon to be enriched dataset
data_target = rem_ND(data_target, 7976)
data_codes = data_target[,1] ##Keeping the genome codes aside for further use


data_target = Enrich(data_target, data_codes, AllSubsets) #Create all new individuals

#Particular case : Colistin, 2 unique individuals