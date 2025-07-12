#!/usr/bin/env python3
"""
Research Paper Search Tool
Integrates with multiple paper databases and AI-powered analysis
"""

import sys
import json
import os
import requests
import time
from datetime import datetime, timedelta
from typing import Dict, List, Optional, Any
import re
from urllib.parse import quote_plus
import xml.etree.ElementTree as ET
from dataclasses import dataclass
import hashlib

# Add current directory to path for imports
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

@dataclass
class Paper:
    """Represents a research paper with metadata"""
    title: str
    authors: List[str]
    abstract: str
    doi: Optional[str] = None
    arxiv_id: Optional[str] = None
    pmid: Optional[str] = None
    url: Optional[str] = None
    published_date: Optional[str] = None
    journal: Optional[str] = None
    keywords: List[str] = None
    score: float = 0.0
    source: str = "unknown"
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert paper to dictionary"""
        return {
            'title': self.title,
            'authors': self.authors,
            'abstract': self.abstract,
            'doi': self.doi,
            'arxiv_id': self.arxiv_id,
            'pmid': self.pmid,
            'url': self.url,
            'published_date': self.published_date,
            'journal': self.journal,
            'keywords': self.keywords or [],
            'score': self.score,
            'source': self.source
        }

class PaperSearchEngine:
    """Comprehensive paper search engine integrating multiple databases"""
    
    def __init__(self):
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'ChemSearch-PaperTool/1.0 (research; contact@example.com)'
        })
        self.cache = {}
        self.rate_limits = {
            'pubmed': 3,  # requests per second
            'arxiv': 3,
            'crossref': 10
        }
        self.last_request = {}
    
    def _rate_limit(self, service: str):
        """Implement rate limiting for API calls"""
        if service in self.last_request:
            elapsed = time.time() - self.last_request[service]
            min_interval = 1.0 / self.rate_limits.get(service, 1)
            if elapsed < min_interval:
                time.sleep(min_interval - elapsed)
        self.last_request[service] = time.time()
    
    def _get_cache_key(self, query: str, source: str, filters: Dict = None) -> str:
        """Generate cache key for search results"""
        cache_data = f"{query}_{source}_{str(filters or {})}"
        return hashlib.md5(cache_data.encode()).hexdigest()
    
    def search_arxiv(self, query: str, max_results: int = 20, category: str = None) -> List[Paper]:
        """Search ArXiv for chemistry and materials science papers"""
        try:
            self._rate_limit('arxiv')
            
            # Build search query
            search_query = query
            if category:
                search_query = f"cat:{category} AND ({query})"
            else:
                # Default to relevant categories for chemistry
                categories = ["cond-mat.*", "physics.chem-ph", "q-bio.*", "cs.LG"]
                cat_query = " OR ".join([f"cat:{cat}" for cat in categories])
                search_query = f"({cat_query}) AND ({query})"
            
            # ArXiv API parameters
            params = {
                'search_query': search_query,
                'start': 0,
                'max_results': max_results,
                'sortBy': 'relevance',
                'sortOrder': 'descending'
            }
            
            url = "http://export.arxiv.org/api/query"
            response = self.session.get(url, params=params, timeout=30)
            response.raise_for_status()
            
            # Parse XML response
            root = ET.fromstring(response.content)
            papers = []
            
            for entry in root.findall('{http://www.w3.org/2005/Atom}entry'):
                try:
                    title = entry.find('{http://www.w3.org/2005/Atom}title').text.strip()
                    abstract = entry.find('{http://www.w3.org/2005/Atom}summary').text.strip()
                    
                    # Extract authors
                    authors = []
                    for author in entry.findall('{http://www.w3.org/2005/Atom}author'):
                        name = author.find('{http://www.w3.org/2005/Atom}name').text
                        authors.append(name)
                    
                    # Extract ArXiv ID and URL
                    arxiv_id = None
                    url = None
                    for link in entry.findall('{http://www.w3.org/2005/Atom}link'):
                        if link.get('title') == 'pdf':
                            url = link.get('href')
                            # Extract arXiv ID from URL
                            match = re.search(r'arxiv\.org/pdf/([^.]+)', url)
                            if match:
                                arxiv_id = match.group(1)
                    
                    # Extract publication date
                    published = entry.find('{http://www.w3.org/2005/Atom}published')
                    published_date = published.text[:10] if published is not None else None
                    
                    # Extract categories as keywords
                    keywords = []
                    for category in entry.findall('{http://arxiv.org/schemas/atom}category', 
                                                namespaces={'': 'http://arxiv.org/schemas/atom'}):
                        keywords.append(category.get('term'))
                    
                    paper = Paper(
                        title=title,
                        authors=authors,
                        abstract=abstract,
                        arxiv_id=arxiv_id,
                        url=url,
                        published_date=published_date,
                        keywords=keywords,
                        source="arxiv"
                    )
                    papers.append(paper)
                    
                except Exception as e:
                    print(f"Error parsing ArXiv entry: {e}", file=sys.stderr)
                    continue
            
            return papers
            
        except Exception as e:
            print(f"Error searching ArXiv: {e}", file=sys.stderr)
            return []
    
    def search_pubmed(self, query: str, max_results: int = 20, retmax: int = None) -> List[Paper]:
        """Search PubMed for biomedical and chemistry papers"""
        try:
            self._rate_limit('pubmed')
            
            if retmax is None:
                retmax = max_results
            
            # Step 1: Search for paper IDs
            search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
            search_params = {
                'db': 'pubmed',
                'term': query,
                'retmax': retmax,
                'sort': 'relevance',
                'retmode': 'json'
            }
            
            response = self.session.get(search_url, params=search_params, timeout=30)
            response.raise_for_status()
            search_data = response.json()
            
            if 'esearchresult' not in search_data or 'idlist' not in search_data['esearchresult']:
                return []
            
            pmids = search_data['esearchresult']['idlist']
            if not pmids:
                return []
            
            # Step 2: Fetch detailed information
            self._rate_limit('pubmed')
            
            fetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
            fetch_params = {
                'db': 'pubmed',
                'id': ','.join(pmids),
                'retmode': 'xml'
            }
            
            response = self.session.get(fetch_url, params=fetch_params, timeout=30)
            response.raise_for_status()
            
            # Parse XML response
            root = ET.fromstring(response.content)
            papers = []
            
            for article in root.findall('.//PubmedArticle'):
                try:
                    # Extract title
                    title_elem = article.find('.//ArticleTitle')
                    title = title_elem.text if title_elem is not None else "Unknown Title"
                    
                    # Extract abstract
                    abstract_elem = article.find('.//Abstract/AbstractText')
                    abstract = abstract_elem.text if abstract_elem is not None else "No abstract available"
                    
                    # Extract authors
                    authors = []
                    for author in article.findall('.//Author'):
                        last_name = author.find('.//LastName')
                        first_name = author.find('.//ForeName')
                        if last_name is not None and first_name is not None:
                            authors.append(f"{first_name.text} {last_name.text}")
                    
                    # Extract PMID
                    pmid_elem = article.find('.//PMID')
                    pmid = pmid_elem.text if pmid_elem is not None else None
                    
                    # Extract DOI
                    doi = None
                    for article_id in article.findall('.//ArticleId'):
                        if article_id.get('IdType') == 'doi':
                            doi = article_id.text
                            break
                    
                    # Extract journal
                    journal_elem = article.find('.//Journal/Title')
                    journal = journal_elem.text if journal_elem is not None else None
                    
                    # Extract publication date
                    pub_date = None
                    date_elem = article.find('.//PubDate')
                    if date_elem is not None:
                        year = date_elem.find('.//Year')
                        month = date_elem.find('.//Month')
                        day = date_elem.find('.//Day')
                        if year is not None:
                            pub_date = year.text
                            if month is not None:
                                pub_date += f"-{month.text.zfill(2)}"
                                if day is not None:
                                    pub_date += f"-{day.text.zfill(2)}"
                    
                    # Extract keywords/MeSH terms
                    keywords = []
                    for mesh in article.findall('.//MeshHeading/DescriptorName'):
                        keywords.append(mesh.text)
                    
                    # Create URL
                    url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else None
                    
                    paper = Paper(
                        title=title,
                        authors=authors,
                        abstract=abstract,
                        pmid=pmid,
                        doi=doi,
                        url=url,
                        published_date=pub_date,
                        journal=journal,
                        keywords=keywords,
                        source="pubmed"
                    )
                    papers.append(paper)
                    
                except Exception as e:
                    print(f"Error parsing PubMed entry: {e}", file=sys.stderr)
                    continue
            
            return papers
            
        except Exception as e:
            print(f"Error searching PubMed: {e}", file=sys.stderr)
            return []
    
    def search_crossref(self, query: str, max_results: int = 20) -> List[Paper]:
        """Search CrossRef for academic papers with DOIs"""
        try:
            self._rate_limit('crossref')
            
            url = "https://api.crossref.org/works"
            params = {
                'query': query,
                'rows': max_results,
                'sort': 'relevance',
                'filter': 'type:journal-article'
            }
            
            response = self.session.get(url, params=params, timeout=30)
            response.raise_for_status()
            data = response.json()
            
            papers = []
            
            for item in data.get('message', {}).get('items', []):
                try:
                    title = item.get('title', ['Unknown Title'])[0]
                    
                    # Extract authors
                    authors = []
                    for author in item.get('author', []):
                        given = author.get('given', '')
                        family = author.get('family', '')
                        if given and family:
                            authors.append(f"{given} {family}")
                    
                    # Extract abstract (often not available in CrossRef)
                    abstract = item.get('abstract', 'Abstract not available from CrossRef')
                    
                    # Extract DOI
                    doi = item.get('DOI')
                    
                    # Extract journal
                    journal = item.get('container-title', [''])[0]
                    
                    # Extract publication date
                    pub_date = None
                    date_parts = item.get('published-print', {}).get('date-parts') or \
                                item.get('published-online', {}).get('date-parts')
                    if date_parts and date_parts[0]:
                        date_list = date_parts[0]
                        if len(date_list) >= 1:
                            pub_date = str(date_list[0])
                            if len(date_list) >= 2:
                                pub_date += f"-{str(date_list[1]).zfill(2)}"
                                if len(date_list) >= 3:
                                    pub_date += f"-{str(date_list[2]).zfill(2)}"
                    
                    # Extract subject areas as keywords
                    keywords = item.get('subject', [])
                    
                    # Create URL
                    url = f"https://doi.org/{doi}" if doi else item.get('URL')
                    
                    paper = Paper(
                        title=title,
                        authors=authors,
                        abstract=abstract,
                        doi=doi,
                        url=url,
                        published_date=pub_date,
                        journal=journal,
                        keywords=keywords,
                        source="crossref"
                    )
                    papers.append(paper)
                    
                except Exception as e:
                    print(f"Error parsing CrossRef entry: {e}", file=sys.stderr)
                    continue
            
            return papers
            
        except Exception as e:
            print(f"Error searching CrossRef: {e}", file=sys.stderr)
            return []
    
    def search_multiple_sources(self, query: str, sources: List[str] = None, 
                              max_results_per_source: int = 10) -> Dict[str, List[Paper]]:
        """Search multiple paper databases simultaneously"""
        if sources is None:
            sources = ['arxiv', 'pubmed', 'crossref']
        
        results = {}
        
        for source in sources:
            print(f"Searching {source}...", file=sys.stderr)
            try:
                if source == 'arxiv':
                    papers = self.search_arxiv(query, max_results_per_source)
                elif source == 'pubmed':
                    papers = self.search_pubmed(query, max_results_per_source)
                elif source == 'crossref':
                    papers = self.search_crossref(query, max_results_per_source)
                else:
                    papers = []
                
                results[source] = papers
                print(f"Found {len(papers)} papers from {source}", file=sys.stderr)
                
            except Exception as e:
                print(f"Error searching {source}: {e}", file=sys.stderr)
                results[source] = []
        
        return results
    
    def deduplicate_papers(self, papers: List[Paper]) -> List[Paper]:
        """Remove duplicate papers based on title and DOI similarity"""
        unique_papers = []
        seen_titles = set()
        seen_dois = set()
        
        for paper in papers:
            # Normalize title for comparison
            normalized_title = re.sub(r'\W+', ' ', paper.title.lower()).strip()
            
            # Check for duplicates
            is_duplicate = False
            
            # Check DOI
            if paper.doi and paper.doi in seen_dois:
                is_duplicate = True
            
            # Check title similarity
            if not is_duplicate:
                for seen_title in seen_titles:
                    # Simple similarity check (can be improved with fuzzy matching)
                    if len(normalized_title) > 10 and normalized_title in seen_title:
                        is_duplicate = True
                        break
            
            if not is_duplicate:
                unique_papers.append(paper)
                seen_titles.add(normalized_title)
                if paper.doi:
                    seen_dois.add(paper.doi)
        
        return unique_papers
    
    def rank_papers(self, papers: List[Paper], query: str) -> List[Paper]:
        """Rank papers by relevance to query"""
        query_terms = set(query.lower().split())
        
        for paper in papers:
            score = 0.0
            
            # Title relevance (highest weight)
            title_terms = set(paper.title.lower().split())
            title_overlap = len(query_terms.intersection(title_terms))
            score += title_overlap * 3.0
            
            # Abstract relevance
            abstract_terms = set(paper.abstract.lower().split())
            abstract_overlap = len(query_terms.intersection(abstract_terms))
            score += abstract_overlap * 1.0
            
            # Keywords relevance
            if paper.keywords:
                keyword_terms = set(' '.join(paper.keywords).lower().split())
                keyword_overlap = len(query_terms.intersection(keyword_terms))
                score += keyword_overlap * 2.0
            
            # Recency bonus (papers from last 3 years get bonus)
            if paper.published_date:
                try:
                    pub_year = int(paper.published_date[:4])
                    current_year = datetime.now().year
                    if current_year - pub_year <= 3:
                        score += 1.0
                except:
                    pass
            
            paper.score = score
        
        return sorted(papers, key=lambda p: p.score, reverse=True)
    
    def search(self, query: str, sources: List[str] = None, max_results: int = 50,
               include_analysis: bool = True) -> Dict[str, Any]:
        """Main search function that combines all sources and provides analysis"""
        try:
            # Search multiple sources
            source_results = self.search_multiple_sources(
                query, sources, max_results // len(sources or ['arxiv', 'pubmed', 'crossref'])
            )
            
            # Combine all papers
            all_papers = []
            for papers in source_results.values():
                all_papers.extend(papers)
            
            # Deduplicate and rank
            unique_papers = self.deduplicate_papers(all_papers)
            ranked_papers = self.rank_papers(unique_papers, query)
            
            # Limit results
            final_papers = ranked_papers[:max_results]
            
            # Prepare results
            results = {
                'success': True,
                'query': query,
                'total_found': len(final_papers),
                'sources_searched': list(source_results.keys()),
                'papers': [paper.to_dict() for paper in final_papers],
                'source_breakdown': {
                    source: len(papers) for source, papers in source_results.items()
                }
            }
            
            # Add analysis if requested
            if include_analysis and final_papers:
                results['analysis'] = self.analyze_papers(final_papers, query)
            
            return results
            
        except Exception as e:
            return {
                'success': False,
                'error': str(e),
                'query': query,
                'papers': []
            }
    
    def analyze_papers(self, papers: List[Paper], query: str) -> Dict[str, Any]:
        """Analyze a collection of papers to extract insights"""
        try:
            analysis = {
                'summary': {
                    'total_papers': len(papers),
                    'date_range': None,
                    'top_journals': [],
                    'top_authors': [],
                    'common_keywords': []
                },
                'trends': {
                    'publication_timeline': {},
                    'research_areas': {}
                },
                'insights': []
            }
            
            # Date range analysis
            dates = [p.published_date for p in papers if p.published_date]
            if dates:
                analysis['summary']['date_range'] = {
                    'earliest': min(dates),
                    'latest': max(dates)
                }
            
            # Journal analysis
            journal_counts = {}
            for paper in papers:
                if paper.journal:
                    journal_counts[paper.journal] = journal_counts.get(paper.journal, 0) + 1
            
            analysis['summary']['top_journals'] = sorted(
                journal_counts.items(), key=lambda x: x[1], reverse=True
            )[:10]
            
            # Author analysis
            author_counts = {}
            for paper in papers:
                for author in paper.authors:
                    author_counts[author] = author_counts.get(author, 0) + 1
            
            analysis['summary']['top_authors'] = sorted(
                author_counts.items(), key=lambda x: x[1], reverse=True
            )[:10]
            
            # Keyword analysis
            keyword_counts = {}
            for paper in papers:
                if paper.keywords:
                    for keyword in paper.keywords:
                        keyword_counts[keyword] = keyword_counts.get(keyword, 0) + 1
            
            analysis['summary']['common_keywords'] = sorted(
                keyword_counts.items(), key=lambda x: x[1], reverse=True
            )[:20]
            
            # Publication timeline
            year_counts = {}
            for paper in papers:
                if paper.published_date:
                    try:
                        year = paper.published_date[:4]
                        year_counts[year] = year_counts.get(year, 0) + 1
                    except:
                        pass
            
            analysis['trends']['publication_timeline'] = dict(sorted(year_counts.items()))
            
            # Generate insights
            insights = []
            
            if len(papers) > 0:
                insights.append(f"Found {len(papers)} relevant papers matching your query.")
            
            if analysis['summary']['date_range']:
                date_range = analysis['summary']['date_range']
                insights.append(f"Papers span from {date_range['earliest'][:4]} to {date_range['latest'][:4]}.")
            
            if analysis['summary']['top_journals']:
                top_journal = analysis['summary']['top_journals'][0]
                insights.append(f"Most papers published in '{top_journal[0]}' ({top_journal[1]} papers).")
            
            if year_counts:
                recent_years = sorted(year_counts.items(), reverse=True)[:3]
                if recent_years[0][1] > 1:
                    insights.append(f"Peak publication activity in {recent_years[0][0]} with {recent_years[0][1]} papers.")
            
            analysis['insights'] = insights
            
            return analysis
            
        except Exception as e:
            print(f"Error in paper analysis: {e}", file=sys.stderr)
            return {'error': str(e)}

def main():
    """Main function for command-line usage"""
    if len(sys.argv) < 2:
        print("Usage: python paper_search.py <query> [sources] [max_results]")
        print("Sources: arxiv,pubmed,crossref (default: all)")
        print("Example: python paper_search.py 'perovskite solar cells' arxiv,pubmed 20")
        sys.exit(1)
    
    query = sys.argv[1]
    sources = sys.argv[2].split(',') if len(sys.argv) > 2 else None
    max_results = int(sys.argv[3]) if len(sys.argv) > 3 else 30
    
    engine = PaperSearchEngine()
    results = engine.search(query, sources, max_results)
    
    # Output JSON results
    print(json.dumps(results, indent=2, ensure_ascii=False))

if __name__ == "__main__":
    main()
