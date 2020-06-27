library(tidyverse) # metapackage with lots of helpful functions
library(purrr)
library(jsonlite)
library(dplyr)
library(tidytext)
library(tidyr)
library(ggplot2)
library(textmineR)
library(magrittr)
library(data.table)
library(topicmodels)
library(colorspace)
library(ldatuning)
library(gmp)
library(wordcloud)
library(RColorBrewer)
library(lubridate)
library(reshape2)
library(tm)
library(writexl)


path <- "C:/Users/qwerty/Downloads/CORD-19-research-challenge/biorxiv_medrxiv/biorxiv_medrxiv"

temp <- list.files(path, pattern="*.json", full.names=TRUE)

df <- purrr::map_df(temp, function(x) { 
  purrr::map(jsonlite::fromJSON(x), function(y) ifelse(is.null(y), NA, y)) 
})


length(df$metadata)
head(df,n=1)
row_index <- rownames(df)
# i will use paper_id and abstract for topic modeling
new_df <- df %>% select(abstract)
new_df$row_id <- rownames(df)
head(new_df,n=5)
view(new_df,n=5)

# Pre-processing-----------------------------

new_df$abstract <- gsub("c\\(","",new_df$abstract)
new_df <- new_df[!(is.na(new_df$row_id) | new_df$row_id=="NULL" | is.na(new_df$abstract) | new_df$abstract=="NULL"),] 

text_cleaning_tokens <- new_df %>% tidytext::unnest_tokens(word, abstract)
words <- text_cleaning_tokens$word
text_cleaning_tokens$word <- gsub('[[:digit:]]+', '', text_cleaning_tokens$word)
text_cleaning_tokens$word <- gsub('[[:punct:]]+', '', text_cleaning_tokens$word)
text_cleaning_tokens <- text_cleaning_tokens %>% filter(!(nchar(word) == 1))%>% anti_join(stop_words)

myStopWords <- c("covid","doi","https","biorxiv","doi_biorxiv","doiorg","https_doiorg","copyright_holder",
               "doiorg_doi","holder","biorxiv_preprint","preprint_peer","copyright","holder_preprint",
               "peer_reviewed","author","author_funder","rights_reserved","reserved_reuse","allowed_permission",
              "permission","funder_rights", "funder","data", "license","virus","viruses","pandemic",
              "license_display","medrxiv_license","preprint_perpetuity", "abstract","word",
              "count","reserved","preprint","text","rights")
text_cleaning_tokens <- text_cleaning_tokens %>% filter(!(text_cleaning_tokens$word %in% myStopWords))
text_cleaning_tokens %>% count(word, sort = TRUE)
tokens <- text_cleaning_tokens %>% filter(!(word==""))
#mutate: create new column called ind
tokens <- tokens %>% mutate(ind = row_number())
#spread() turns a pair of key:value columns
tokens <- tokens %>% group_by(row_id) %>% mutate(ind = row_number()) %>%
  tidyr::spread(key = ind, value = word)
tokens [is.na(tokens)] <- ""
tokens <- tidyr::unite(tokens, abstract,-row_id,sep =" " )
tokens$abstract <- trimws(tokens$abstract)

#create DocumentTermMatrix
dtm <- CreateDtm(tokens$abstract, 
                 doc_names = tokens$row_id, 
                 ngram_window = c(1, 2))

#dtm <- dtm[, !(colnames(dtm) %in% myStopWords)]
#sel_idx <- rowSums(dtm) > 0
#dtm <- dtm[sel_idx, ]
#textdata <- textdata[sel_idx, ]
#explore the basic frequency
tf <- TermDocFreq(dtm = dtm)
original_tf <- tf %>% select(term, term_freq,doc_freq)
rownames(original_tf) <- 1:nrow(original_tf)
# Eliminate words appearing less than 2 times or in more than half of the
# documents
vocabulary <- tf$term[ tf$term_freq > 1 & tf$doc_freq < nrow(dtm) / 2 ]
dtm <- dtm[ , vocabulary]
#dim(dtm)
dtm = dtm

# Running LDA -----------------------------------------------------------
k_list <- seq(2, 30, by = 1)
model_dir <- paste0("models_", digest::digest(vocabulary, algo = "sha1"))
if (!dir.exists(model_dir)) dir.create(model_dir)
model_list <- TmParallelApply(X = k_list, FUN = function(k){
  filename = file.path(model_dir, paste0(k, "_topics.rda"))
  
  if (!file.exists(filename)) {
    model_lda <- FitLdaModel(dtm = dtm, k = k, iterations = 500)
    model_lda$k <- k
    model_lda$coherence <- CalcProbCoherence(phi = model_lda$phi, dtm = dtm, M = 5)
    save(model_lda, file = filename)
  } else {
    load(filename)
  }
  
  model_lda
}, export=c("dtm", "model_dir"))
#model tuning
#choosing the best model
coherence_mat <- data.frame(k = sapply(model_list, function(x) nrow(x$phi)), 
                            coherence = sapply(model_list, function(x) mean(x$coherence)), 
                            stringsAsFactors = FALSE)
ggplot(coherence_mat, aes(x = k, y = coherence)) +
  geom_point() +
  geom_line(group = 1)+
  ggtitle("Best Topic by Coherence Score") + theme_minimal() +
  scale_x_continuous(breaks = seq(2,30,1)) + ylab("Coherence")

plot(coherence_mat, type = "o")

# select k based on maximum average coherence
best_model <- model_list[which.max(coherence_mat$coherence)][[ 1 ]]
#model <- model_list[ coherence_mat$coherence == max(coherence_mat$coherence) ][[ 1 ]]
names(best_model) # phi is P(words | topics), theta is P(topics | documents)

# Get the R-squared of this model
best_model$r2 <- CalcTopicModelR2(dtm = dtm, phi = best_model$phi, theta = best_model$theta)
best_model$r2

# probabilistic coherence, a measure of topic quality
# this measure can be used with any topic model, not just probabilistic ones
summary(best_model$coherence)

hist(best_model$coherence,
     col="blue",
     main = "Histogram of probabilistic coherence")

#Top 20 terms based on phi value (pr(word|topic)) — probability of word given a topic 
best_model$top_terms <- GetTopTerms(phi = best_model$phi, M = 20)
top20_wide <- as.data.frame(best_model$top_terms)
write_xlsx(top20_wide,"top20_new6.xlsx")
head(t(best_model$top_terms))
# Probabilistic coherence: measures statistical support for a topic
best_model$coherence <- CalcProbCoherence(phi = best_model$phi, dtm = dtm, M = 5)
best_model_coherence <- as.data.frame(best_model$coherence)
write_xlsx(best_model_coherence,"coherence_6.xlsx")

# Get the prevalence of each topic
best_model$prevalence <- colSums(best_model$theta) / sum(best_model$theta) * 100

# prevalence should be proportional to alpha
plot(best_model$prevalence, best_model$alpha, xlab = "prevalence", ylab = "alpha")


# textmineR has a naive topic labeling tool based on probable bigrams
best_model$labels <- LabelTopics(assignments = best_model$theta > 0.05, 
                            dtm = dtm,
                            M = 1)

head(best_model$labels)
#summary of 10 most prevalent topics
best_model$summary <- data.frame(topic = rownames(best_model$phi),
                                 label = best_model$labels,
                                 coherence = round(best_model$coherence,3),
                                 prevalence = round(best_model$prevalence,3),
                                 top_terms = apply(best_model$top_terms,2, function(x){
                                   paste(x,collapse = ",")
                                 }),
                                 stringsAsFactors = FALSE)

prevalent_topics <- best_model$summary[order(best_model$summary$prevalence, decreasing = TRUE), ][1:10, ]
write_xlsx(prevalent_topics,"prevalent_topics_6.xlsx")

#prediction with gibbs
predict_gibbs <- predict(best_model, dtm,
                         method = "gibbs",
                         iterations = 200,
                         burnin = 180,
                         cpus = 2)

#predictions with dot
predict_dot <- predict(best_model, dtm,
                       method = "dot")

#comparison
barplot(rbind(predict_gibbs[10,],predict_dot[10,]),
        col = c("red","blue"), las = 2, beside = TRUE)
legend("topleft",legend = c("gibbs","dot"), col = c("red","blue"),
       fill = c("red","blue"),cex = 0.5)


# word, topic relationship  pr(word|topic) ---------------------------------------------

allterms <-data.frame(t(best_model$phi))
allterms$word <- rownames(allterms)
rownames(allterms) <- 1:nrow(allterms)
allterms <- melt(allterms,idvars = "word") 
allterms <- allterms %>% rename(topic = variable)
FINAL_allterms <- allterms %>% group_by(topic) %>% arrange(desc(value))

#Topic,word,freq ------------------------------------------------------
final_summary_words <- data.frame(top_terms = t(best_model$top_terms))
final_summary_words$topic <- rownames(final_summary_words)
rownames(final_summary_words) <- 1:nrow(final_summary_words)
final_summary_words <- final_summary_words %>% melt(id.vars = c("topic"))
final_summary_words <- final_summary_words %>% rename(word = value) %>% select(-variable)
final_summary_words <- left_join(final_summary_words,allterms)
final_summary_words <- final_summary_words %>% group_by(topic,word) %>%
  arrange(desc(value))
final_summary_words <- final_summary_words %>% group_by(topic, word) %>% filter(row_number() == 1) %>% 
  ungroup() %>% tidyr::separate(topic, into =c("t","topic")) %>% select(-t)
word_topic_freq <- left_join(final_summary_words, original_tf, by = c("word" = "term"))
write_xlsx(word_topic_freq,"word_topic_freq_6.xlsx")

#per-document-per-topic probabilities ----------------------------------------------

theta_df <- data.frame(best_model$theta)
theta_df$document <-rownames(theta_df) 
rownames(theta_df) <- 1:nrow(theta_df)
theta_df$document <- as.numeric(theta_df$document)
theta_df <- melt(theta_df,id.vars = "document")
theta_df <- theta_df %>% rename(topic = variable) 
theta_df <- theta_df %>% tidyr::separate(topic, into =c("t","topic")) %>% select(-t)
FINAL_document_topic <- theta_df %>% group_by(document) %>% 
  arrange(desc(value)) %>% filter(row_number() ==1)
write_xlsx(FINAL_document_topic,"FINAL_document_topic_6.xlsx")

#Visualising of topics in a dendrogram (Hellinger distance) ----------------------------------------------

best_model$topic_linguistic_dist <- CalcHellingerDist(best_model$phi)
best_model$hclust <- hclust(as.dist(best_model$topic_linguistic_dist), "ward.D")
best_model$hclust$labels <- paste(best_model$hclust$labels, best_model$labels[ , 1])
plot(best_model$hclust)

#visualising topics of words based on the max value of phi
set.seed(1234)
pdf("cluster6.pdf")
for(i in 1:length(unique(final_summary_words$topic)))

  {  wordcloud(words = subset(final_summary_words ,topic == i)$word, freq = subset(final_summary_words ,topic == i)$value, min.freq = 1,
             max.words=200, random.order=FALSE, rot.per=0.35, 
             colors=brewer.pal(8, "Dark2"))}
dev.off()

################### LSA MODEL ##################

# get a tf-idf matrix
tf_sample <- TermDocFreq(dtm)

tf_sample$idf[ is.infinite(tf_sample$idf) ] <- 0 # fix idf for missing words

tf_idf <- t(dtm / rowSums(dtm)) * tf_sample$idf

tf_idf <- t(tf_idf)

# Fit a Latent Semantic Analysis model
lsa_model <- FitLsaModel(dtm = tf_idf, 
                         k = 30)

# probabilistic coherence, a measure of topic quality

summary(lsa_model$coherence)

hist(lsa_model$coherence, col= "blue")

# Get the top terms of each topic
lsa_model$top_terms <- GetTopTerms(phi = lsa_model$phi, M = 5)
head(t(lsa_model$top_terms))
lsa_top_terms <- as.data.frame(t(lsa_model$top_terms))
write_xlsx(lsa_top_terms,"lsa_model_top_terms.xlsx")

# Get the prevalence of each topic
lsa_model$prevalence <- colSums(lsa_model$theta) / sum(lsa_model$theta) * 100

# textmineR has a naive topic labeling tool based on probable bigrams
lsa_model$labels <- LabelTopics(assignments = lsa_model$theta > 0.05, 
                                dtm = dtm,
                                M = 5)
head(lsa_model$labels)

# coherence into a summary table
lsa_model$summary <- data.frame(topic = rownames(lsa_model$phi),
                                label = lsa_model$labels,
                                coherence = round(lsa_model$coherence, 3),
                                prevalence = round(lsa_model$prevalence,3),
                                top_terms = apply(lsa_model$top_terms, 2, function(x){
                                  paste(x, collapse = ", ")
                                }),
                                stringsAsFactors = FALSE)
lsa_prevalent_topics <- lsa_model$summary[ order(lsa_model$summary$prevalence, decreasing = TRUE) , ][ 1:10 , ]

write_xlsx(lsa_prevalent_topics,"lsa_prevalent_topics.xlsx")

#  predictions
lsa_predict <- t(dtm / rowSums(dtm)) * tf_sample$idf

lsa_predict <- t(lsa_predict)

lsa_predict <- predict(lsa_model, lsa_predict)

# comparison
barplot(rbind(lsa_model$theta[ rownames(dtm)[ 1 ] , ],
              lsa_predict[ rownames(dtm)[ 1 ] , ]), 
        las = 2,
        main = "Comparing topic predictions in LSA",
        beside = TRUE,
        col = c("red", "blue"))

legend("topright", 
       legend = c("During fitting", "Predicted"),
       fill = c("red", "blue"), cex = 0.5)

lsa_gamma <- as.data.frame(CalcGamma(phi=lsa_model$phi,
                   theta = lsa_model$theta))
lsa_model$top_terms <- GetTopTerms(phi = lsa_model$phi, M = 10)
top10_lsa <- as.data.frame(lsa_model$top_terms)
write_xlsx(top10_lsa,"top10_lsa.xlsx")

top10_lsa_transpose <- as.data.frame(t(top10_lsa))
top10_lsa_transpose$topicsterms <- paste(top10_lsa_transpose$V1, top10_lsa_transpose$V2,
                               top10_lsa_transpose$V3, top10_lsa_transpose$V4,
                               top10_lsa_transpose$V5,top10_lsa_transpose$V6,
                               top10_lsa_transpose$V7, top10_lsa_transpose$V8,
                               top10_lsa_transpose$V9, top10_lsa_transpose$V10)

top10_lsa_transpose$topiclabel <- rownames(top10_lsa_transpose)

topics_lsa_df <- top10_lsa_transpose[,-(1:10),drop = FALSE]

rownames(topics_lsa_df) <- NULL

####################   RECOMMENDATION SYSTEM ##########

####### LDA RECOMMENDATION ######
top20_wide_transpose <- as.data.frame(t(top20_wide))
top20_wide_transpose$topicsterms <- paste(top20_wide_transpose$V1,top20_wide_transpose$V2,
                                          top20_wide_transpose$V3, top20_wide_transpose$V4,
                                          top20_wide_transpose$V5,top20_wide_transpose$V6,
                                          top20_wide_transpose$V7, top20_wide_transpose$V8,
                                          top20_wide_transpose$V9, top20_wide_transpose$V10)
top20_wide_transpose$topiclabel <- rownames(top20_wide_transpose)

topics_df <- top20_wide_transpose[,-(1:20),drop = FALSE]

rownames(topics_df) <- NULL

topics_df$rowIndex <- as.numeric(rownames(topics_df))

topicList <- as.list(topics_df$topicsterms)

N.topics <- length(topicList)

QuerySearch <- function(queryTerm) {
  
# Record starting time 
start.time <- Sys.time()

# store docs in Corpus class 
my.topics <- VectorSource(c(topicList, queryTerm))

# Transform/standaridze topics 
my.corpus <- VCorpus(my.topics) %>% 
              tm_map(stemDocument) %>%
              tm_map(removeNumbers) %>% 
              tm_map(content_transformer(tolower)) %>% 
              tm_map(removeWords,stopwords("en")) %>%
              tm_map(stripWhitespace)

# Store docs into a term document matrix 
term.topics.matrix.stm <- TermDocumentMatrix(my.corpus,
                                         control=list(
                                         weighting=function(x) weightSMART(x,spec="ltc"),
                                         wordLengths=c(1,Inf)))

# Transform term document matrix into a dataframe
term.topics.matrix <- tidy(term.topics.matrix.stm) %>% 
                      group_by(document) %>% 
                      mutate(vtrLen=sqrt(sum(count^2))) %>% 
                      mutate(count=count/vtrLen) %>% 
                      ungroup() %>% 
                      select(term:count)

topicsMatrix <- term.topics.matrix %>% 
                mutate(document=as.numeric(document)) %>% 
                filter(document<N.topics+1)

qryMatrix <- term.topics.matrix %>% 
              mutate(document=as.numeric(document)) %>% 
              filter(document>=N.topics+1)

# Calcualte results by cosine similarity
searchResult <- topicsMatrix %>% 
              inner_join(qryMatrix,by=c("term"="term"),
              suffix=c(".doc",".query")) %>% 
              mutate(termScore=round(count.doc*count.query,4)) %>% 
              group_by(document.query,document.doc) %>% 
              summarise(Score=sum(termScore)) %>% 
              filter(row_number(desc(Score))<=10) %>% 
              arrange(desc(Score)) %>% 
              left_join(topics_df,by=c("document.doc"="rowIndex")) %>% 
              ungroup() %>% 
              rename(Result=topicsterms) %>% 
              select(Result,Score,topiclabel) %>% 
              data.frame()

# Record when it stops and take the difference
end.time <- Sys.time()
time.taken <- round(end.time - start.time,4)
print(paste("Used",time.taken,"seconds"))

return(searchResult)

}

############## LSA recommendatýon ####

top10_lsa_transpose <- as.data.frame(t(top10_lsa))
top10_lsa_transpose$topicsterms <- paste(top10_lsa_transpose$V1, top10_lsa_transpose$V2,
                                         top10_lsa_transpose$V3, top10_lsa_transpose$V4,
                                         top10_lsa_transpose$V5,top10_lsa_transpose$V6,
                                         top10_lsa_transpose$V7, top10_lsa_transpose$V8,
                                         top10_lsa_transpose$V9, top10_lsa_transpose$V10)

top10_lsa_transpose$topiclabel <- rownames(top10_lsa_transpose)

topics_lsa_df <- top10_lsa_transpose[,-(1:10),drop = FALSE]

rownames(topics_lsa_df) <- NULL

topics_lsa_df$rowIndex <- as.numeric(rownames(topics_lsa_df))

topicList_lsa <- as.list(topics_lsa_df$topicsterms)

N.topics_lsa <- length(topicList_lsa)

QuerySearch_LSA <- function(queryTerm) {
  
  # Record starting time 
  start.time <- Sys.time()
  
  # store docs in Corpus class
  my.topics_lsa <- VectorSource(c(topicList_lsa, queryTerm))
  
  # Transform/standaridze topics 
  my.corpus_lsa <- VCorpus(my.topics_lsa) %>% 
    tm_map(stemDocument) %>%
    tm_map(removeNumbers) %>% 
    tm_map(content_transformer(tolower)) %>% 
    tm_map(removeWords,stopwords("en")) %>%
    tm_map(stripWhitespace)
  
  # Store docs into a term document matrix 
  term.topics.matrix.stm_lsa <- TermDocumentMatrix(my.corpus_lsa,
                                               control=list(
                                                 weighting=function(x) weightSMART(x,spec="ltc"),
                                                 wordLengths=c(1,Inf)))
  
  # Transform term document matrix into a dataframe
  term.topics.matrix_lsa <- tidy(term.topics.matrix.stm_lsa) %>% 
    group_by(document) %>% 
    mutate(vtrLen=sqrt(sum(count^2))) %>% 
    mutate(count=count/vtrLen) %>% 
    ungroup() %>% 
    select(term:count)
  
  topicsMatrix_lsa <- term.topics.matrix_lsa %>% 
    mutate(document=as.numeric(document)) %>% 
    filter(document<N.topics_lsa + 1)
  
  qryMatrix <- term.topics.matrix_lsa %>% 
    mutate(document=as.numeric(document)) %>% 
    filter(document>=N.topics_lsa+1)
  
  # Calcualte top results by cosine similarity
  searchResult_lsa <- topicsMatrix_lsa %>% 
    inner_join(qryMatrix,by=c("term"="term"),
               suffix=c(".doc",".query")) %>% 
    mutate(termScore=round(count.doc*count.query,4)) %>% 
    group_by(document.query,document.doc) %>% 
    summarise(Score=sum(termScore)) %>% 
    filter(row_number(desc(Score))<=10) %>% 
    arrange(desc(Score)) %>% 
    left_join(topics_df,by=c("document.doc"="rowIndex")) %>% 
    ungroup() %>% 
    rename(Result=topicsterms) %>% 
    select(Result,Score,topiclabel) %>% 
    data.frame()
  
  # Record when it stops and take the difference
  end.time <- Sys.time()
  time.taken <- round(end.time - start.time,4)
  print(paste("Used",time.taken,"seconds"))
  
  
  return(searchResult_lsa)
  
}

QuerySearch_LSA ("patients clinical china covid")
QuerySearch("patients clinical china covid")
