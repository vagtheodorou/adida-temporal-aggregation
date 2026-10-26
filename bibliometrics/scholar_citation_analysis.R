################################################################################
# scholar_citation_analysis.R
# Bibliometric analysis of the works citing the original ADIDA paper
# (Nikolopoulos et al., 2011): Google Scholar "Cited by" scraping, abstract
# enrichment via OpenAlex / Crossref / Semantic Scholar, lemmatised corpus,
# word frequencies and word cloud (Section 2 of the paper).
################################################################################

# Run from the folder where the outputs should be written.

# Required packages
required_packages <- c(
  "scholar", "rvest", "httr", "jsonlite", "tm", "wordcloud",
  "RColorBrewer", "dplyr", "stringr", "tidyr", "xml2", "textstem"
)

new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

lapply(required_packages, library, character.only = TRUE)

################################################################################
# PART 1: EXTRACT CITATIONS WITH CITATION COUNTS
################################################################################

#' Extract citing papers from Google Scholar with citation counts
#' 
#' NEW: Now extracts how many times each citing paper has been cited
extract_scholar_citations_with_counts <- function(cite_id, max_pages = 10) {
  
  cat("Extracting citations from Google Scholar...\n")
  cat("NEW: Also extracting citation counts for each paper\n\n")
  
  all_citations <- list()
  
  for (page in 0:(max_pages - 1)) {
    
    start <- page * 10
    url <- paste0(
      "https://scholar.google.com/scholar?",
      "cites=", cite_id,
      "&as_sdt=2005&sciodt=0,5&hl=en",
      "&start=", start
    )
    
    cat(sprintf("Scraping page %d (start=%d)...\n", page + 1, start))
    
    if (page > 0) Sys.sleep(runif(1, 2, 4))
    
    tryCatch({
      page_html <- read_html(url)
      results <- page_html %>% html_nodes(".gs_ri")
      
      if (length(results) == 0) {
        cat("No more results found. Stopping.\n")
        break
      }
      
      for (result in results) {
        
        # Title
        title_node <- result %>% html_node(".gs_rt")
        title <- title_node %>% html_text() %>% str_trim()
        title <- gsub("\\[PDF\\]|\\[HTML\\]|\\[BOOK\\]", "", title) %>% str_trim()
        
        # Link
        link <- title_node %>% html_node("a") %>% html_attr("href")
        
        # Authors and venue info
        info_node <- result %>% html_node(".gs_a")
        info_text <- info_node %>% html_text()
        info_parts <- str_split(info_text, " - ")[[1]]
        
        authors <- if(length(info_parts) >= 1) info_parts[1] else NA
        year_match <- str_extract(info_text, "\\d{4}")
        year <- if(!is.na(year_match)) as.integer(year_match) else NA
        venue <- if(length(info_parts) >= 2) info_parts[2] else NA
        
        # Snippet
        snippet <- result %>% html_node(".gs_rs") %>% html_text() %>% str_trim()
        
        # NEW: Extract citation count
        # Look for "Cited by X" link in .gs_fl (footer links)
        cited_by <- result %>% html_nodes(".gs_fl a:nth-child(3)") %>% html_text()
        cited_by <- as.numeric(str_extract(cited_by, "\\d+"))
        
        all_citations[[length(all_citations) + 1]] <- data.frame(
          title = title,
          authors = authors,
          year = year,
          venue = venue,
          link = ifelse(is.na(link), NA, link),
          snippet = snippet,
          cited_by = cited_by,  # NEW FIELD
          stringsAsFactors = FALSE
        )
      }
      
    }, error = function(e) {
      cat(sprintf("Error on page %d: %s\n", page + 1, e$message))
    })
  }
  
  if (length(all_citations) == 0) {
    stop("No citations extracted. Check the cite_id or Google Scholar availability.")
  }
  
  citations_df <- bind_rows(all_citations)
  
  # Summary statistics
  cat(sprintf("\nExtracted %d citing papers from Google Scholar\n", nrow(citations_df)))
  cat(sprintf("Citation count range: %d to %d\n", 
              min(citations_df$cited_by, na.rm = TRUE),
              max(citations_df$cited_by, na.rm = TRUE)))
  cat(sprintf("Mean citations per paper: %.1f\n", 
              mean(citations_df$cited_by, na.rm = TRUE)))
  
  # Show top cited papers
  cat("\nTop 5 most cited papers:\n")
  top_cited <- citations_df %>% 
    arrange(desc(cited_by)) %>% 
    head(5) %>%
    select(title, cited_by)
  print(top_cited)
  cat("\n")
  
  return(citations_df)
}


################################################################################
# PART 2: IMPROVED API ENRICHMENT (with httr:: namespace)
################################################################################

query_openalex_improved <- function(title, debug = TRUE) {
  
  if (is.na(title) || nchar(title) < 5) return(NULL)
  
  clean_title <- title
  clean_title <- gsub("[[:punct:]]", " ", clean_title)
  clean_title <- str_squish(clean_title)
  
  search_query <- URLencode(clean_title)
  url <- paste0("https://api.openalex.org/works?search=", search_query)
  
  if (debug) cat("  -> Querying OpenAlex...\n")
  
  tryCatch({
    response <- GET(url, 
                    add_headers(`User-Agent` = "mailto:researcher@example.com"),
                    timeout(10))
    
    if (status_code(response) != 200) {
      if (debug) cat("     OpenAlex returned status:", status_code(response), "\n")
      return(NULL)
    }
    
    data <- httr::content(response, "parsed", encoding = "UTF-8")
    
    if (is.null(data$results) || length(data$results) == 0) {
      if (debug) cat("     No results from OpenAlex\n")
      return(NULL)
    }
    
    work <- data$results[[1]]
    result_title <- work$title
    
    if (!is.null(result_title)) {
      title_lower <- tolower(clean_title)
      result_lower <- tolower(result_title)
      
      if (!grepl(substr(title_lower, 1, 20), result_lower) && 
          !grepl(substr(result_lower, 1, 20), title_lower)) {
        if (debug) cat("     Title mismatch, skipping\n")
        return(NULL)
      }
    }
    
    abstract_text <- NA
    if (!is.null(work$abstract_inverted_index) && length(work$abstract_inverted_index) > 0) {
      abstract_text <- reconstruct_abstract(work$abstract_inverted_index)
      if (debug && !is.na(abstract_text)) {
        cat("     ✓ Found abstract (", nchar(abstract_text), " chars)\n", sep="")
      }
    }
    
    doi <- if (!is.null(work$doi)) work$doi else NA
    oa_url <- NULL
    is_oa <- FALSE
    if (!is.null(work$open_access)) {
      oa_url <- work$open_access$oa_url
      is_oa <- isTRUE(work$open_access$is_oa)
    }
    
    return(list(
      doi = doi,
      abstract = abstract_text,
      title = result_title,
      oa_url = oa_url,
      is_oa = is_oa,
      source = "openalex"
    ))
    
  }, error = function(e) {
    if (debug) cat("     OpenAlex error:", e$message, "\n")
    return(NULL)
  })
}

reconstruct_abstract <- function(inverted_index) {
  if (is.null(inverted_index) || length(inverted_index) == 0) return(NA)
  
  tryCatch({
    pairs <- list()
    for (word in names(inverted_index)) {
      positions <- inverted_index[[word]]
      for (pos in positions) {
        pairs[[as.character(pos)]] <- word
      }
    }
    
    if (length(pairs) == 0) return(NA)
    
    positions <- as.integer(names(pairs))
    sorted_indices <- order(positions)
    sorted_words <- unlist(pairs)[sorted_indices]
    
    abstract <- paste(sorted_words, collapse = " ")
    abstract <- str_squish(abstract)
    
    if (nchar(abstract) < 50) return(NA)
    
    return(abstract)
  }, error = function(e) {
    return(NA)
  })
}

query_crossref_improved <- function(title, debug = TRUE) {
  
  if (is.na(title) || nchar(title) < 5) return(NULL)
  
  clean_title <- gsub("[[:punct:]]", " ", title)
  clean_title <- str_squish(clean_title)
  
  search_query <- URLencode(clean_title)
  url <- paste0("https://api.crossref.org/works?query.title=", search_query, "&rows=1")
  
  if (debug) cat("  -> Querying Crossref...\n")
  
  tryCatch({
    response <- GET(url, timeout(10))
    
    if (status_code(response) != 200) {
      if (debug) cat("     Crossref returned status:", status_code(response), "\n")
      return(NULL)
    }
    
    data <- httr::content(response, "parsed", encoding = "UTF-8")
    
    if (is.null(data$message$items) || length(data$message$items) == 0) {
      if (debug) cat("     No results from Crossref\n")
      return(NULL)
    }
    
    item <- data$message$items[[1]]
    
    abstract_text <- NA
    if (!is.null(item$abstract)) {
      abstract_clean <- gsub("<.*?>", "", item$abstract)
      abstract_clean <- gsub("&lt;.*?&gt;", "", abstract_clean)
      abstract_clean <- str_squish(abstract_clean)
      
      if (nchar(abstract_clean) > 50) {
        abstract_text <- abstract_clean
        if (debug) cat("     ✓ Found abstract (", nchar(abstract_text), " chars)\n", sep="")
      }
    } else {
      if (debug) cat("     No abstract in Crossref record\n")
    }
    
    doi <- if (!is.null(item$DOI)) item$DOI else NA
    result_title <- if (!is.null(item$title)) paste(item$title, collapse = " ") else NA
    
    return(list(
      doi = doi,
      abstract = abstract_text,
      title = result_title,
      source = "crossref"
    ))
    
  }, error = function(e) {
    if (debug) cat("     Crossref error:", e$message, "\n")
    return(NULL)
  })
}

query_semantic_scholar <- function(title, debug = TRUE) {
  
  if (is.na(title) || nchar(title) < 5) return(NULL)
  
  clean_title <- gsub("[[:punct:]]", " ", title)
  clean_title <- str_squish(clean_title)
  
  search_query <- URLencode(clean_title)
  url <- paste0("https://api.semanticscholar.org/graph/v1/paper/search?query=", 
                search_query, "&fields=title,abstract,openAccessPdf")
  
  if (debug) cat("  -> Querying Semantic Scholar...\n")
  
  tryCatch({
    response <- GET(url, timeout(10))
    
    if (status_code(response) != 200) {
      if (debug) cat("     Semantic Scholar returned status:", status_code(response), "\n")
      return(NULL)
    }
    
    data <- httr::content(response, "parsed", encoding = "UTF-8")
    
    if (is.null(data$data) || length(data$data) == 0) {
      if (debug) cat("     No results from Semantic Scholar\n")
      return(NULL)
    }
    
    paper <- data$data[[1]]
    
    abstract_text <- NA
    if (!is.null(paper$abstract) && nchar(paper$abstract) > 50) {
      abstract_text <- paper$abstract
      if (debug) cat("     ✓ Found abstract (", nchar(abstract_text), " chars)\n", sep="")
    } else {
      if (debug) cat("     No abstract in Semantic Scholar record\n")
    }
    
    return(list(
      abstract = abstract_text,
      title = paper$title,
      source = "semantic_scholar"
    ))
    
  }, error = function(e) {
    if (debug) cat("     Semantic Scholar error:", e$message, "\n")
    return(NULL)
  })
}

enrich_citations_improved <- function(citations_df, debug = TRUE) {
  
  cat("\n=== IMPROVED ENRICHMENT PROCESS ===\n")
  cat("Trying multiple APIs: OpenAlex → Crossref → Semantic Scholar\n")
  cat("This may take several minutes...\n\n")
  
  citations_df$abstract <- NA_character_
  citations_df$full_text <- NA_character_
  citations_df$text_source <- "none"
  citations_df$doi <- NA_character_
  
  for (i in 1:nrow(citations_df)) {
    
    cat(sprintf("\n[%d/%d] %s (cited by %d)\n", 
                i, nrow(citations_df), 
                substr(citations_df$title[i], 1, 60),
                citations_df$cited_by[i]))
    
    found_abstract <- FALSE
    
    Sys.sleep(0.3)
    openalex_data <- query_openalex_improved(citations_df$title[i], debug = debug)
    
    if (!is.null(openalex_data)) {
      if (!is.null(openalex_data$doi)) {
        citations_df$doi[i] <- openalex_data$doi
      }
      
      if (!is.null(openalex_data$abstract) && !is.na(openalex_data$abstract)) {
        if (nchar(openalex_data$abstract) > 50) {
          citations_df$abstract[i] <- openalex_data$abstract
          citations_df$text_source[i] <- "abstract_openalex"
          found_abstract <- TRUE
          cat("  ✓✓ SUCCESS: Got abstract from OpenAlex\n")
        }
      }
    }
    
    if (!found_abstract) {
      Sys.sleep(0.3)
      crossref_data <- query_crossref_improved(citations_df$title[i], debug = debug)
      
      if (!is.null(crossref_data)) {
        if (!is.null(crossref_data$doi) && is.na(citations_df$doi[i])) {
          citations_df$doi[i] <- crossref_data$doi
        }
        
        if (!is.null(crossref_data$abstract) && !is.na(crossref_data$abstract)) {
          if (nchar(crossref_data$abstract) > 50) {
            citations_df$abstract[i] <- crossref_data$abstract
            citations_df$text_source[i] <- "abstract_crossref"
            found_abstract <- TRUE
            cat("  ✓✓ SUCCESS: Got abstract from Crossref\n")
          }
        }
      }
    }
    
    if (!found_abstract) {
      Sys.sleep(0.3)
      semantic_data <- query_semantic_scholar(citations_df$title[i], debug = debug)
      
      if (!is.null(semantic_data)) {
        if (!is.null(semantic_data$abstract) && !is.na(semantic_data$abstract)) {
          if (nchar(semantic_data$abstract) > 50) {
            citations_df$abstract[i] <- semantic_data$abstract
            citations_df$text_source[i] <- "abstract_semantic"
            found_abstract <- TRUE
            cat("  ✓✓ SUCCESS: Got abstract from Semantic Scholar\n")
          }
        }
      }
    }
    
    if (!found_abstract && !is.na(citations_df$snippet[i])) {
      if (nchar(citations_df$snippet[i]) > 20) {
        citations_df$abstract[i] <- citations_df$snippet[i]
        citations_df$text_source[i] <- "snippet_scholar"
        cat("  → Using Google Scholar snippet (no abstract found)\n")
      }
    }
    
    if (!found_abstract && citations_df$text_source[i] == "none") {
      cat("  ✗ No text found - will use title only\n")
    }
  }
  
  cat("\n\n=== ENRICHMENT SUMMARY ===\n")
  cat(sprintf("Total papers: %d\n", nrow(citations_df)))
  cat(sprintf("With abstracts (OpenAlex): %d (%.1f%%)\n", 
              sum(citations_df$text_source == "abstract_openalex"),
              100 * sum(citations_df$text_source == "abstract_openalex") / nrow(citations_df)))
  cat(sprintf("With abstracts (Crossref): %d (%.1f%%)\n", 
              sum(citations_df$text_source == "abstract_crossref"),
              100 * sum(citations_df$text_source == "abstract_crossref") / nrow(citations_df)))
  cat(sprintf("With abstracts (Semantic Scholar): %d (%.1f%%)\n", 
              sum(citations_df$text_source == "abstract_semantic"),
              100 * sum(citations_df$text_source == "abstract_semantic") / nrow(citations_df)))
  cat(sprintf("With snippets only: %d (%.1f%%)\n", 
              sum(citations_df$text_source == "snippet_scholar"),
              100 * sum(citations_df$text_source == "snippet_scholar") / nrow(citations_df)))
  cat(sprintf("No text: %d (%.1f%%)\n", 
              sum(citations_df$text_source == "none"),
              100 * sum(citations_df$text_source == "none") / nrow(citations_df)))
  
  total_with_abstracts <- sum(citations_df$text_source %in% 
                                c("abstract_openalex", "abstract_crossref", "abstract_semantic"))
  cat(sprintf("\n✓ TOTAL WITH FULL ABSTRACTS: %d/%d (%.1f%%)\n", 
              total_with_abstracts, nrow(citations_df),
              100 * total_with_abstracts / nrow(citations_df)))
  
  return(citations_df)
}


################################################################################
# PART 3: TEXT CLEANING (with lemmatization)
################################################################################

create_corpus <- function(citations_df) {
  
  cat("\n\nCreating text corpus...\n")
  
  corpus_texts <- character(nrow(citations_df))
  
  for (i in 1:nrow(citations_df)) {
    text_parts <- c()
    
    if (!is.na(citations_df$title[i])) {
      text_parts <- c(text_parts, citations_df$title[i])
    }
    
    if (!is.na(citations_df$abstract[i])) {
      text_parts <- c(text_parts, citations_df$abstract[i])
    }
    
    corpus_texts[i] <- paste(text_parts, collapse = " ")
  }
  
  corpus_texts <- corpus_texts[nchar(corpus_texts) > 0]
  
  cat(sprintf("Corpus created with %d documents\n", length(corpus_texts)))
  
  return(corpus_texts)
}

clean_corpus_improved <- function(corpus_texts) {
  
  cat("\n=== IMPROVED TEXT CLEANING ===\n")
  cat("Using lemmatization (not stemming) for better word forms\n")
  
  corpus <- VCorpus(VectorSource(corpus_texts))
  
  cat("  - Converting to lowercase...\n")
  corpus <- tm_map(corpus, content_transformer(tolower))
  
  cat("  - Removing punctuation...\n")
  corpus <- tm_map(corpus, removePunctuation)
  
  cat("  - Removing numbers...\n")
  corpus <- tm_map(corpus, removeNumbers)
  
  cat("  - Removing standard stopwords...\n")
  corpus <- tm_map(corpus, removeWords, stopwords("english"))
  
  cat("  - Removing academic boilerplate...\n")
  academic_stopwords <- c(
    "study", "studies", "research", "paper", "article", "articles",
    "results", "result", "method", "methods", "methodology", "data", "analysis",
    "used", "using", "use", "uses", "show", "shows", "showed", "shown",
    "approach", "approaches", "model", "models", "modeling", "modelling",
    "based", "propose", "proposed", "proposes", "however", "also",
    "can", "may", "will", "shall", "one", "two", "three", "first", "second",
    "new", "different", "important", "significant", "significantly",
    "table", "figure", "et", "al", "elsevier", "rights", "reserved",
    "university", "press", "journal", "publisher", "doi", "httpsdoi",
    "concluded", "conclude", "found", "find", "findings",
    "presented", "present", "investigated", "investigate",
    "demonstrated", "demonstrate", "examined", "examine"
  )
  
  corpus <- tm_map(corpus, removeWords, academic_stopwords)
  
  cat("  - Removing extra whitespace...\n")
  corpus <- tm_map(corpus, stripWhitespace)
  
  cat("  - Applying lemmatization (this preserves readable words)...\n")
  corpus <- tm_map(corpus, content_transformer(function(x) {
    lemmatize_strings(x)
  }))
  
  cat("✓ Text cleaning complete!\n")
  
  return(corpus)
}


################################################################################
# PART 4: WORD FREQUENCY AND WORD CLOUD
################################################################################

get_word_frequencies <- function(corpus) {
  
  cat("\nCalculating word frequencies...\n")
  
  dtm <- DocumentTermMatrix(corpus)
  freq <- colSums(as.matrix(dtm))
  freq <- sort(freq, decreasing = TRUE)
  
  word_freq <- data.frame(
    word = names(freq),
    freq = freq,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  
  cat(sprintf("Found %d unique terms\n", nrow(word_freq)))
  cat("\nTop 20 most frequent words:\n")
  print(head(word_freq, 20))
  
  return(word_freq)
}

generate_wordcloud <- function(word_freq, output_file = "citation_wordcloud.png", 
                               max_words = 100) {
  
  cat(sprintf("\nGenerating word cloud with top %d words...\n", max_words))
  
  png(output_file, width = 1200, height = 800, res = 150)
  
  colors <- brewer.pal(8, "Dark2")
  
  wordcloud(
    words = word_freq$word,
    freq = word_freq$freq,
    max.words = max_words,
    random.order = FALSE,
    rot.per = 0.35,
    colors = colors,
    scale = c(4, 0.5)
  )
  
  title(main = "Word Cloud: Papers Citing My Work",
        sub = paste("Based on", nrow(word_freq), "unique terms from citing papers"),
        cex.main = 1.5)
  
  dev.off()
  
  cat(sprintf("Word cloud saved to: %s\n", output_file))
}


################################################################################
# MAIN PIPELINE
################################################################################

main_pipeline_final <- function(cite_id = "15938694505943547654", 
                                max_pages = 10,
                                output_file = "citation_wordcloud_final.png",
                                debug = TRUE,
                                min_citations = 0) {  # NEW: Filter by citation count
  
  cat("================================================================================\n")
  cat("FINAL GOOGLE SCHOLAR CITATION WORD CLOUD PIPELINE\n")
  cat("================================================================================\n\n")
  
  cat("FEATURES:\n")
  cat("  ✓ Extracts citation counts for each paper\n")
  cat("  ✓ Multiple API sources (OpenAlex + Crossref + Semantic Scholar)\n")
  cat("  ✓ Lemmatization for readable words\n")
  cat("  ✓ Optional filtering by citation count\n\n")
  
  # Step 1: Extract citations WITH COUNTS
  citations <- extract_scholar_citations_with_counts(cite_id, max_pages)
  write.csv(citations, "citations_raw.csv", row.names = FALSE)
  cat("\nRaw citations saved to: citations_raw.csv\n")
  
  # Optional: Filter by citation count
  if (min_citations > 0) {
    cat(sprintf("\nFiltering papers with at least %d citations...\n", min_citations))
    n_before <- nrow(citations)
    citations <- citations %>% filter(cited_by >= min_citations)
    cat(sprintf("Kept %d out of %d papers (%.1f%%)\n", 
                nrow(citations), n_before,
                100 * nrow(citations) / n_before))
  }
  
  # Step 2: Enrichment
  citations_enriched <- enrich_citations_improved(citations, debug = debug)
  write.csv(citations_enriched, "citations_enriched.csv", row.names = FALSE)
  cat("\nEnriched citations saved to: citations_enriched.csv\n")
  
  # Step 3: Create and clean corpus
  corpus_texts <- create_corpus(citations_enriched)
  corpus_clean <- clean_corpus_improved(corpus_texts)
  
  # Step 4: Word frequencies
  word_freq <- get_word_frequencies(corpus_clean)
  write.csv(word_freq, "word_frequencies.csv", row.names = FALSE)
  cat("\nWord frequencies saved to: word_frequencies.csv\n")
  
  # Step 5: Word cloud
  generate_wordcloud(word_freq, output_file)
  
  cat("\n================================================================================\n")
  cat("PIPELINE COMPLETE!\n")
  cat("================================================================================\n")
  cat("\nOutputs:\n")
  cat("  1. citations_raw.csv - WITH CITATION COUNTS\n")
  cat("  2. citations_enriched.csv - WITH ABSTRACTS + CITATION COUNTS\n")
  cat("  3. word_frequencies.csv - READABLE WORDS\n")
  cat("  4.", output_file, "- Word cloud\n\n")
  
  cat("ANALYSIS IDEAS:\n")
  cat("  - Sort by cited_by to find most influential papers\n")
  cat("  - Filter by min_citations for focused analysis\n")
  cat("  - Compare word clouds of highly-cited vs less-cited papers\n\n")
  
  invisible(list(
    citations = citations_enriched,
    word_freq = word_freq
  ))
}


################################################################################
# RUN THE PIPELINE
################################################################################

# BASIC RUN - All papers
results <- main_pipeline_final(
  cite_id = "15938694505943547654",
  max_pages = 30,
  output_file = "citation_wordcloud_final.png",
  debug = TRUE
)

# OPTIONAL: View citations sorted by influence
cat("\n\n=== VIEW CITATIONS SORTED BY INFLUENCE ===\n")
citations_sorted <- read.csv("citations_enriched.csv") %>%
  arrange(desc(cited_by)) %>%
  select(title, authors, year, cited_by, text_source)

cat("\nTop 10 most cited papers:\n")
print(head(citations_sorted, 10))

cat("\n✓✓✓ COMPLETE! Check citations_raw.csv for citation counts ✓✓✓\n")


### WORDCLOUD EDITING ###

word_freq<-read.csv("word_frequencies.csv")
max_words = 200
colors <- brewer.pal(8, "Paired")
wordcloud(
  words = word_freq$word,
  freq = word_freq$freq,
  max.words = max_words,
  random.order = FALSE,
  rot.per = 0.35,
  colors = colors,
  scale = c(4, 0.5)
)

title(main = "Word Cloud: Papers Citing My Work",
      sub = paste("Based on", nrow(word_freq), "unique terms from citing papers"),
      cex.main = 1.5)