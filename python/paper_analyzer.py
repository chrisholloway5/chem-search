#!/usr/bin/env python3
"""
AI-Powered Paper Analysis Tool
Provides advanced analysis capabilities for research papers
"""

import sys
import json
import os
from typing import Dict, List, Any, Optional
import re
from collections import defaultdict, Counter
import nltk
from nltk.corpus import stopwords
from nltk.tokenize import word_tokenize, sent_tokenize
from nltk.stem import WordNetLemmatizer
import openai
from datetime import datetime

# Download required NLTK data (only if not already present)
try:
    nltk.data.find('tokenizers/punkt')
except LookupError:
    nltk.download('punkt', quiet=True)

try:
    nltk.data.find('corpora/stopwords')
except LookupError:
    nltk.download('stopwords', quiet=True)

try:
    nltk.data.find('corpora/wordnet')
except LookupError:
    nltk.download('wordnet', quiet=True)

# Add current directory to path for imports
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

class PaperAnalyzer:
    """Advanced AI-powered paper analysis tool"""
    
    def __init__(self):
        self.stop_words = set(stopwords.words('english'))
        self.lemmatizer = WordNetLemmatizer()
        
        # Chemistry-specific terms and patterns
        self.chemistry_terms = {
            'synthesis': ['synthesis', 'synthesize', 'synthetic', 'preparation', 'fabrication'],
            'catalysis': ['catalyst', 'catalysis', 'catalytic', 'catalyze', 'reaction'],
            'materials': ['material', 'polymer', 'composite', 'nanoparticle', 'crystal'],
            'characterization': ['characterization', 'analysis', 'spectroscopy', 'microscopy'],
            'properties': ['property', 'properties', 'behavior', 'performance', 'stability'],
            'drug_discovery': ['drug', 'pharmaceutical', 'therapeutic', 'medicine', 'compound'],
            'energy': ['battery', 'solar', 'fuel', 'energy', 'electrochemical', 'photovoltaic'],
            'environment': ['environmental', 'green', 'sustainable', 'pollution', 'remediation']
        }
        
        # Load OpenAI API key if available
        self.openai_client = None
        try:
            openai_key = os.getenv('OPENAI_API_KEY')
            if openai_key:
                openai.api_key = openai_key
                self.openai_client = openai
        except:
            pass
    
    def extract_key_findings(self, papers: List[Dict]) -> List[Dict[str, Any]]:
        """Extract key findings from paper abstracts"""
        findings = []
        
        for paper in papers:
            try:
                abstract = paper.get('abstract', '')
                if not abstract or len(abstract) < 50:
                    continue
                
                # Extract sentences that likely contain findings
                sentences = sent_tokenize(abstract)
                finding_sentences = []
                
                # Look for sentences with finding indicators
                finding_indicators = [
                    'show', 'demonstrate', 'reveal', 'indicate', 'suggest',
                    'find', 'observe', 'report', 'achieve', 'obtain',
                    'result', 'conclude', 'prove', 'confirm'
                ]
                
                for sentence in sentences:
                    sentence_lower = sentence.lower()
                    if any(indicator in sentence_lower for indicator in finding_indicators):
                        finding_sentences.append(sentence)
                
                if finding_sentences:
                    finding = {
                        'paper_title': paper.get('title', 'Unknown'),
                        'authors': paper.get('authors', []),
                        'key_sentences': finding_sentences[:3],  # Top 3 findings
                        'source': paper.get('source', 'unknown'),
                        'score': paper.get('score', 0)
                    }
                    findings.append(finding)
                    
            except Exception as e:
                print(f"Error extracting findings from paper: {e}", file=sys.stderr)
                continue
        
        # Sort by paper score
        return sorted(findings, key=lambda x: x['score'], reverse=True)
    
    def identify_research_methods(self, papers: List[Dict]) -> Dict[str, List[str]]:
        """Identify common research methods across papers"""
        methods = defaultdict(list)
        
        # Method indicators
        method_patterns = {
            'synthesis_methods': [
                r'hydrothermal', r'sol-gel', r'chemical vapor deposition', r'cvd',
                r'electrodeposition', r'solvothermal', r'microwave', r'sonochemical'
            ],
            'characterization_methods': [
                r'x-ray diffraction', r'xrd', r'nmr', r'ftir', r'raman',
                r'sem', r'tem', r'afm', r'xps', r'uv-vis'
            ],
            'computational_methods': [
                r'dft', r'density functional theory', r'molecular dynamics',
                r'monte carlo', r'machine learning', r'neural network'
            ],
            'analytical_methods': [
                r'hplc', r'gc-ms', r'mass spectrometry', r'chromatography',
                r'electrochemical', r'cyclic voltammetry'
            ]
        }
        
        for paper in papers:
            text = (paper.get('title', '') + ' ' + paper.get('abstract', '')).lower()
            
            for method_type, patterns in method_patterns.items():
                for pattern in patterns:
                    if re.search(pattern, text):
                        methods[method_type].append(paper.get('title', 'Unknown'))
        
        # Remove duplicates and count occurrences
        for method_type in methods:
            methods[method_type] = list(set(methods[method_type]))
        
        return dict(methods)
    
    def analyze_research_trends(self, papers: List[Dict]) -> Dict[str, Any]:
        """Analyze research trends over time"""
        trends = {
            'temporal_trends': {},
            'emerging_topics': [],
            'declining_topics': [],
            'hot_topics': []
        }
        
        try:
            # Group papers by year
            papers_by_year = defaultdict(list)
            for paper in papers:
                pub_date = paper.get('published_date')
                if pub_date:
                    try:
                        year = int(pub_date[:4])
                        papers_by_year[year].append(paper)
                    except:
                        continue
            
            # Analyze terms by year
            terms_by_year = defaultdict(lambda: defaultdict(int))
            
            for year, year_papers in papers_by_year.items():
                for paper in year_papers:
                    text = (paper.get('title', '') + ' ' + paper.get('abstract', '')).lower()
                    words = word_tokenize(text)
                    words = [self.lemmatizer.lemmatize(word) for word in words 
                            if word.isalpha() and word not in self.stop_words and len(word) > 3]
                    
                    for word in words:
                        terms_by_year[year][word] += 1
            
            # Identify trending terms
            if len(papers_by_year) >= 2:
                years = sorted(papers_by_year.keys())
                recent_years = years[-2:]  # Last 2 years
                older_years = years[:-2] if len(years) > 2 else []
                
                if older_years and recent_years:
                    # Compare recent vs older periods
                    recent_terms = defaultdict(int)
                    older_terms = defaultdict(int)
                    
                    for year in recent_years:
                        for term, count in terms_by_year[year].items():
                            recent_terms[term] += count
                    
                    for year in older_years:
                        for term, count in terms_by_year[year].items():
                            older_terms[term] += count
                    
                    # Find emerging topics (high in recent, low in older)
                    emerging = []
                    for term, recent_count in recent_terms.items():
                        older_count = older_terms.get(term, 0)
                        if recent_count >= 3 and (older_count == 0 or recent_count / older_count > 2):
                            emerging.append((term, recent_count))
                    
                    trends['emerging_topics'] = sorted(emerging, key=lambda x: x[1], reverse=True)[:10]
            
            # Identify hot topics (frequently mentioned overall)
            all_terms = defaultdict(int)
            for paper in papers:
                text = (paper.get('title', '') + ' ' + paper.get('abstract', '')).lower()
                words = word_tokenize(text)
                words = [self.lemmatizer.lemmatize(word) for word in words 
                        if word.isalpha() and word not in self.stop_words and len(word) > 3]
                
                for word in words:
                    all_terms[word] += 1
            
            trends['hot_topics'] = sorted(all_terms.items(), key=lambda x: x[1], reverse=True)[:20]
            
            return trends
            
        except Exception as e:
            print(f"Error analyzing trends: {e}", file=sys.stderr)
            return trends
    
    def classify_research_areas(self, papers: List[Dict]) -> Dict[str, List[str]]:
        """Classify papers into research areas"""
        classifications = defaultdict(list)
        
        for paper in papers:
            title = paper.get('title', '').lower()
            abstract = paper.get('abstract', '').lower()
            text = title + ' ' + abstract
            
            # Check against chemistry term categories
            for area, terms in self.chemistry_terms.items():
                if any(term in text for term in terms):
                    classifications[area].append(paper.get('title', 'Unknown'))
        
        # Remove duplicates
        for area in classifications:
            classifications[area] = list(set(classifications[area]))
        
        return dict(classifications)
    
    def generate_research_gaps(self, papers: List[Dict], query: str) -> List[str]:
        """Identify potential research gaps"""
        gaps = []
        
        try:
            # Analyze what's missing based on common research patterns
            methods = self.identify_research_methods(papers)
            areas = self.classify_research_areas(papers)
            
            # Check for underrepresented areas
            if len(papers) > 5:
                if 'computational_methods' in methods and len(methods['computational_methods']) < len(papers) * 0.3:
                    gaps.append("Computational modeling and simulations could be more widely applied")
                
                if 'characterization_methods' in methods and len(methods['characterization_methods']) < len(papers) * 0.4:
                    gaps.append("More comprehensive characterization techniques could provide better insights")
                
                # Check for missing application areas
                if 'drug_discovery' not in areas and 'pharmaceutical' in query.lower():
                    gaps.append("Limited focus on pharmaceutical applications")
                
                if 'environment' not in areas and any(term in query.lower() for term in ['green', 'sustainable', 'environmental']):
                    gaps.append("Environmental applications and sustainability aspects need more attention")
            
            # Generic research gaps
            gaps.extend([
                "Long-term stability studies and practical implementation challenges",
                "Cost-effectiveness analysis and scalability considerations",
                "Interdisciplinary approaches combining multiple research areas"
            ])
            
        except Exception as e:
            print(f"Error generating research gaps: {e}", file=sys.stderr)
        
        return gaps[:5]  # Return top 5 gaps
    
    def generate_ai_summary(self, papers: List[Dict], query: str) -> Optional[str]:
        """Generate AI-powered summary using OpenAI"""
        if not self.openai_client or len(papers) == 0:
            return None
        
        try:
            # Prepare context from top papers
            context_papers = papers[:5]  # Use top 5 papers
            context = f"Research Query: {query}\n\n"
            
            for i, paper in enumerate(context_papers, 1):
                context += f"Paper {i}:\n"
                context += f"Title: {paper.get('title', 'Unknown')}\n"
                context += f"Abstract: {paper.get('abstract', 'No abstract')[:500]}...\n\n"
            
            prompt = f"""Based on the following research papers related to "{query}", provide a comprehensive summary that includes:

1. Main research themes and approaches
2. Key findings and breakthroughs
3. Current challenges and limitations
4. Future research directions
5. Practical applications and implications

{context}

Please provide a well-structured summary in 300-400 words."""
            
            response = self.openai_client.ChatCompletion.create(
                model="gpt-3.5-turbo",
                messages=[
                    {"role": "system", "content": "You are a research analyst specializing in chemistry and materials science. Provide clear, accurate, and insightful analysis of research papers."},
                    {"role": "user", "content": prompt}
                ],
                max_tokens=500,
                temperature=0.7
            )
            
            return response.choices[0].message.content.strip()
            
        except Exception as e:
            print(f"Error generating AI summary: {e}", file=sys.stderr)
            return None
    
    def analyze_citation_potential(self, papers: List[Dict]) -> List[Dict[str, Any]]:
        """Analyze papers for citation potential and impact"""
        analysis = []
        
        for paper in papers:
            try:
                score = 0
                factors = []
                
                # Check publication date (recency bonus)
                pub_date = paper.get('published_date')
                if pub_date:
                    try:
                        pub_year = int(pub_date[:4])
                        current_year = datetime.now().year
                        years_ago = current_year - pub_year
                        if years_ago <= 1:
                            score += 3
                            factors.append("Very recent publication")
                        elif years_ago <= 3:
                            score += 2
                            factors.append("Recent publication")
                    except:
                        pass
                
                # Check for methodology papers
                title = paper.get('title', '').lower()
                abstract = paper.get('abstract', '').lower()
                
                if any(term in title + abstract for term in ['method', 'technique', 'approach', 'protocol']):
                    score += 2
                    factors.append("Methodological contribution")
                
                # Check for review characteristics
                if any(term in title for term in ['review', 'survey', 'overview', 'perspective']):
                    score += 3
                    factors.append("Review/survey paper")
                
                # Check for novelty indicators
                if any(term in title + abstract for term in ['novel', 'new', 'first', 'breakthrough']):
                    score += 2
                    factors.append("Novelty claims")
                
                # Check journal quality (simplified)
                journal = paper.get('journal', '').lower()
                high_impact_keywords = ['nature', 'science', 'cell', 'journal', 'proceedings']
                if any(keyword in journal for keyword in high_impact_keywords):
                    score += 2
                    factors.append("High-impact journal")
                
                analysis.append({
                    'title': paper.get('title', 'Unknown'),
                    'citation_potential_score': score,
                    'factors': factors,
                    'recommendation': 'High potential' if score >= 5 else 'Medium potential' if score >= 3 else 'Standard'
                })
                
            except Exception as e:
                print(f"Error analyzing citation potential: {e}", file=sys.stderr)
                continue
        
        return sorted(analysis, key=lambda x: x['citation_potential_score'], reverse=True)
    
    def comprehensive_analysis(self, papers: List[Dict], query: str) -> Dict[str, Any]:
        """Perform comprehensive analysis of papers"""
        try:
            analysis = {
                'query': query,
                'total_papers': len(papers),
                'timestamp': datetime.now().isoformat(),
                'key_findings': self.extract_key_findings(papers),
                'research_methods': self.identify_research_methods(papers),
                'research_areas': self.classify_research_areas(papers),
                'trends_analysis': self.analyze_research_trends(papers),
                'research_gaps': self.generate_research_gaps(papers, query),
                'citation_analysis': self.analyze_citation_potential(papers),
                'ai_summary': self.generate_ai_summary(papers, query)
            }
            
            return analysis
            
        except Exception as e:
            return {
                'error': str(e),
                'query': query,
                'timestamp': datetime.now().isoformat()
            }

def main():
    """Main function for command-line usage"""
    if len(sys.argv) < 2:
        print("Usage: python paper_analyzer.py <papers_json> [query]")
        print("Example: python paper_analyzer.py papers.json 'machine learning chemistry'")
        sys.exit(1)
    
    papers_file = sys.argv[1]
    query = sys.argv[2] if len(sys.argv) > 2 else "chemistry research"
    
    try:
        with open(papers_file, 'r', encoding='utf-8') as f:
            data = json.load(f)
        
        papers = data.get('papers', [])
        if not papers:
            print(json.dumps({'error': 'No papers found in input file'}))
            sys.exit(1)
        
        analyzer = PaperAnalyzer()
        analysis = analyzer.comprehensive_analysis(papers, query)
        
        print(json.dumps(analysis, indent=2, ensure_ascii=False))
        
    except Exception as e:
        print(json.dumps({'error': str(e)}))
        sys.exit(1)

if __name__ == "__main__":
    main()
