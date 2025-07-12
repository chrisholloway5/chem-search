#!/usr/bin/env python3
"""
Paper Search Tool Demo
Demonstrates the complete functionality of the research paper search system
"""

import sys
import json
import os
from paper_search import PaperSearchEngine
from paper_analyzer import PaperAnalyzer

def demo_paper_search():
    """Demonstrate the paper search functionality"""
    print("=== Research Paper Search Tool Demo ===\n")
    
    # Initialize search engine
    engine = PaperSearchEngine()
    analyzer = PaperAnalyzer()
    
    # Define test queries for different research areas
    demo_queries = [
        {
            'query': 'perovskite solar cells efficiency',
            'description': 'Energy Research - Solar Cell Materials',
            'sources': ['arxiv', 'pubmed']
        },
        {
            'query': 'machine learning drug discovery QSAR',
            'description': 'Drug Discovery - AI Applications',
            'sources': ['arxiv', 'crossref']
        },
        {
            'query': 'MOF metal organic frameworks CO2 capture',
            'description': 'Environmental Chemistry - Carbon Capture',
            'sources': ['arxiv']
        }
    ]
    
    for i, demo in enumerate(demo_queries, 1):
        print(f"Demo {i}: {demo['description']}")
        print(f"Query: '{demo['query']}'")
        print(f"Sources: {', '.join(demo['sources'])}")
        print("-" * 60)
        
        try:
            # Perform search
            results = engine.search(
                query=demo['query'],
                sources=demo['sources'],
                max_results=10,
                include_analysis=True
            )
            
            if results['success']:
                print(f"✓ Found {results['total_found']} papers")
                
                # Show top papers
                print("\nTop 3 Papers:")
                for j, paper in enumerate(results['papers'][:3], 1):
                    print(f"{j}. {paper['title']}")
                    print(f"   Authors: {', '.join(paper['authors'][:3])}{'...' if len(paper['authors']) > 3 else ''}")
                    print(f"   Source: {paper['source'].upper()}")
                    print(f"   Score: {paper['score']:.1f}")
                    print()
                
                # Perform additional analysis
                if results['papers']:
                    print("Performing AI Analysis...")
                    analysis = analyzer.comprehensive_analysis(results['papers'], demo['query'])
                    
                    if 'research_areas' in analysis and analysis['research_areas']:
                        print("\nResearch Areas Identified:")
                        for area, papers in analysis['research_areas'].items():
                            area_name = area.replace('_', ' ').title()
                            print(f"• {area_name}: {len(papers)} papers")
                    
                    if 'research_gaps' in analysis and analysis['research_gaps']:
                        print(f"\nResearch Gaps Identified:")
                        for gap in analysis['research_gaps'][:3]:
                            print(f"• {gap}")
                    
                    if 'trends_analysis' in analysis and 'hot_topics' in analysis['trends_analysis']:
                        hot_topics = analysis['trends_analysis']['hot_topics'][:5]
                        if hot_topics:
                            print(f"\nHot Topics:")
                            for topic, count in hot_topics:
                                print(f"• {topic} ({count} mentions)")
                
            else:
                print(f"✗ Search failed: {results.get('error', 'Unknown error')}")
        
        except Exception as e:
            print(f"✗ Error: {e}")
        
        print("\n" + "=" * 80 + "\n")
    
    # Summary of capabilities
    print("=== Feature Summary ===")
    print("✓ Multi-database search (ArXiv, PubMed, CrossRef)")
    print("✓ AI-powered relevance ranking")
    print("✓ Automatic deduplication")
    print("✓ Research area classification")
    print("✓ Trend analysis and hot topics")
    print("✓ Research gap identification")
    print("✓ Citation potential analysis")
    print("✓ Key finding extraction")
    print("✓ Method identification")
    print("✓ Export capabilities")
    
    print("\n=== Integration Points ===")
    print("• Web interface: tools/paper_search.php")
    print("• Python backend: python/paper_search.py")
    print("• AI analysis: python/paper_analyzer.py")
    print("• Supports 100+ chemistry research topics")
    print("• Real-time search across academic databases")
    print("• Comprehensive analysis with machine learning")

def demo_specific_features():
    """Demonstrate specific advanced features"""
    print("\n=== Advanced Features Demo ===\n")
    
    engine = PaperSearchEngine()
    
    # Demonstrate different search strategies
    query = "graphene synthesis applications"
    
    print("1. Source-specific search:")
    for source in ['arxiv', 'pubmed', 'crossref']:
        try:
            results = engine.search(query, sources=[source], max_results=5)
            if results['success']:
                print(f"   {source.upper()}: {results['total_found']} papers found")
            else:
                print(f"   {source.upper()}: Search failed")
        except:
            print(f"   {source.upper()}: Error occurred")
    
    print("\n2. Research method identification:")
    try:
        results = engine.search(query, max_results=15)
        if results['success'] and results['papers']:
            analyzer = PaperAnalyzer()
            methods = analyzer.identify_research_methods(results['papers'])
            
            for method_type, papers in methods.items():
                if papers:
                    method_name = method_type.replace('_', ' ').title()
                    print(f"   {method_name}: {len(papers)} papers")
    except Exception as e:
        print(f"   Error in method analysis: {e}")
    
    print("\n3. Citation potential analysis:")
    try:
        if results['success'] and results['papers']:
            citation_analysis = analyzer.analyze_citation_potential(results['papers'])
            print("   Top papers by citation potential:")
            for analysis in citation_analysis[:3]:
                score = analysis['citation_potential_score']
                recommendation = analysis['recommendation']
                title = analysis['title'][:50] + "..." if len(analysis['title']) > 50 else analysis['title']
                print(f"   • {title}")
                print(f"     Score: {score}, Potential: {recommendation}")
    except Exception as e:
        print(f"   Error in citation analysis: {e}")

if __name__ == "__main__":
    # Run demos
    demo_paper_search()
    demo_specific_features()
    
    print("\n" + "=" * 80)
    print("Demo completed! The Research Paper Search Tool is ready for use.")
    print("Access the web interface at: http://localhost/chem-search/tools/paper_search.php")
    print("=" * 80)
