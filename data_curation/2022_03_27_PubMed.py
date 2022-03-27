#!/usr/bin/env python3
"""Query PubMed for related literature"""
import os
from Bio import Entrez, Medline


Entrez.email = "liuxs@shanghaitech.edu.cn"                                      
Entrez.api_key = "6315d840c62b9b7e7928552e51d1e052eb08"


class QueryPubMed:
    def __init__(self, keywords, retmax=1000000):
        self.keywords = keywords
        self.retmax = retmax
        self.count = self.get_count()
        print(self)

    def __repr__(self):
        return f"Search '{self.keywords}', get {self.count} results."

    __str__ = __repr__

    def get_count(self):
        handle = Entrez.egquery(term=self.keywords)
        record = Entrez.read(handle)
        for row in record["eGQueryResult"]:
            if row["DbName"] == "pubmed":
                count = row["Count"]
        return count

    def search(self):                                                           
       handle = Entrez.esearch(db="pubmed", term=self.keywords, retmax=self.retmax)
       record = Entrez.read(handle)                                            
       return record["IdList"]

    
    def save_abs(self, path):
        if not os.path.exists(path):
            os.mkdir(path)
        idlist = self.search()
        idlists = [idlist[i:i+10000] for i in range(0,len(idlist), 10000)]
        i = 0
        count = self.count

        for ids in idlists:
            handle = Entrez.efetch(db="pubmed", id=ids, rettype="medline", retmode="text")
            records = Medline.parse(handle)
            for record in records:
                i += 1
                pmid = record.get("PMID", "")
                print(f"download {pmid}, {i}/{count}")
                article = os.path.join(path, f"{pmid}.txt")
                if os.path.exists(article):
                    continue
                title = record.get("TI", "")
                abstract = record.get("AB", "")
                if not abstract:
                    continue
                with open(article, "w") as f:
                    f.write(f"{title}\n{abstract}\n")

    @staticmethod
    def date(idlist):
        handle = Entrez.efetch(db="pubmed", id=idlist, rettype="medline", retmode="text")
        records = Medline.parse(handle)
        for record in records:
            pmid = record.get("PMID", "")
            date = record.get("DP", "")
            print(pmid, date)



result = QueryPubMed("neopeptide OR cancer peptide vaccine")
result.save_abs("data/")
