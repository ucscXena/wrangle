
#include <map>
#include <iostream>
#include <list>
#include <set>
#include <vector>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

using namespace std;

#define MISSING_VAL -99

class segment {
	public:
	string target;
	string chrome;
	int start, end;
	float value;
};

typedef map<string,list<segment> > chromemap;
typedef map< string, chromemap > segmap;
typedef map<string,set<int> > breakmap;
typedef vector<vector<float> > genemap;


segmap * new_segment() {
	return new segmap();
}

set<string> * new_target_set() {
	return new set<string>();
}

void add_segment(segmap * seg, set<string> *targetSet, char *sample, char *chrom, int chrom_start, int chrom_end, float value) {
	string tmp;
	segment a;
	
	a.target = strdup(sample);
	a.chrome = strdup(chrom);
	a.start = chrom_start;
	a.end = chrom_end;
	a.value = value;
		
	if ( a.target.size() > 0 ) {
		if ( a.chrome.compare("X") == 0 )
			a.chrome = "23";
		if ( a.chrome.compare("Y") == 0 )
			a.chrome = "24";
		
		if ( a.chrome.compare("23") == 0 )
			a.chrome = "chrX";
		else if ( a.chrome.compare("24") == 0 )
			a.chrome = "chrY";
		else if ( a.chrome.find("chr") != 0 )
			a.chrome = string("chr") + a.chrome;
		
		if(seg->find(a.target) == seg->end()) {
			map<string, list<segment> > l;
			(*seg)[ a.target ] = l;
		} 
		if ( (*seg)[ a.target ].find( a.chrome ) == (*seg)[ a.target ].end() ) {
			list<segment> l;
			(*seg)[ a.target ][ a.chrome ] = l;
		}
		(*seg)[ a.target ][a.chrome].push_back( a );
		targetSet->insert( a.target );
	}	
}



void print_matrix(segmap *data, set<string> *targetSet, void (*print)(const char *)) {
	
	//create a break map
	breakmap breaks;
	for( segmap::iterator t = data->begin(); t!=data->end(); ++t) {
		for ( chromemap::iterator c = t->second.begin(); c != t->second.end(); ++c) {
			for ( list<segment>::iterator s = c->second.begin(); s != c->second.end(); s++ ) {
				if ( breaks.find( s->chrome ) == breaks.end() ) {
					set<int> ns;
					ns.insert( s->start );
					ns.insert( s->end );
					breaks[ s->chrome ] = ns;
				} else {
					breaks[ s->chrome ].insert( s->start );
					breaks[ s->chrome ].insert( s->end );
				}
			}
		}	
	}
	
	//print out the matrix column names
	(*print)( "probe" );
	for ( set<string>::iterator ts = targetSet->begin(); ts != targetSet->end(); ++ts ) {
		(*print)( "\t" );
		(*print)( (*ts).c_str() );
	}
	(*print)( "\n" );
	int targetCount = targetSet->size();
	//check breakmap for each chrome
	for ( breakmap::iterator b = breaks.begin(); b != breaks.end(); ++b ) {
		vector<int> starts;
		vector<int> ends;
		vector<string> probeName;
		
		breakmap::iterator cset = breaks.find( b->first );
		int blockCount = cset->second.size() - 1;
		starts.resize( blockCount );
		ends.resize( blockCount );
		probeName.resize( blockCount );
		int i = 0;
		for ( set<int>::iterator v = cset->second.begin(); v != cset->second.end(); ++v ) {
			if ( i > 0 ) {
				ends[i-1] = (*v);
			}
			if ( i < blockCount ) {
				starts[ i ] = (*v);
			}
			i++;
		}
		
		//create names for the breakpoints
		genemap gm;
		gm.resize( blockCount );
		for ( int i = 0; i < blockCount; i++ ) {
			stringstream name;
			name << b->first;
			name << "_";
			name << starts[i];
			name << "_";
			name << ends[i];
			probeName[i] = name.str();
			gm[ i ].resize( targetCount ); 
		}
		//scan the targets		
		int curTarget = 0;
		for ( set<string>::iterator t = targetSet->begin(); t != targetSet->end(); ++t ) {
			//find the target values for the current chrome
			segmap::iterator s = data->find( (*t) );
			chromemap::iterator ss = s->second.find( b->first );
			if ( ss != s->second.end() ) {
				vector<float> vals;
				vals.resize( blockCount );
				for ( int i = 0; i < blockCount; i++ ) {
					vals[ i ] = MISSING_VAL;
				}
				//if the segment overlaps the break, assign the value
				for ( list<segment>::iterator cs = ss->second.begin(); cs != ss->second.end(); cs++ ) {
					for ( int i = 0; i < blockCount; i++ ) {
						if ( cs->end > starts[i] && cs->start < ends[i] ) {
							vals[i] = cs->value;
						}
					}
				}
				//assign the values to the named map
				for ( int i = 0; i < blockCount; i++ ) {
					gm[ i ][ curTarget ] = vals[ i ];
				}				
				//cout << s->first << "\t" << ss->first << "\n";
			}
			curTarget++;
		}
		//print out this chrome's segments
		//and the values for each target
		for ( int i =0; i < blockCount; i++ ) {
			print(probeName[i].c_str());
			for ( int j = 0; j < targetCount; j++ ) {
				float val = gm[i][j];
				print("\t");
				if ( val == MISSING_VAL )
					print("NA");
				else {							
					char str[20]  = "";
					sprintf(str, "%f", val);
					print(str);
				}
			}
			print("\n");
		}
	}
	
}


/*
* This program read's the stdin looking for a segment format of 
* <targetname> <chrome_number> <chrome_start> <chrome_end>
* breaks the segments into discreate probes (so segments on different targets
* have the same probe name) and creates a matrix of probe values for each target
*/


int main(int argc, char** argv) {
	
	segmap data;
	set<string> targetSet;
	
	//load in segment map
	do { 
		string tmp;
		segment a;
		cin >> a.target;
		cin >> a.chrome;
		cin >> tmp;
		a.start = atoi( tmp.c_str() );
		cin >> tmp;
		a.end = atoi( tmp.c_str() );
		cin >> tmp;
		a.value = atof( tmp.c_str() );
		if ( a.target.size() > 0 ) {
			if ( a.chrome.compare("X") == 0 )
				a.chrome = "23";
			if ( a.chrome.compare("Y") == 0 )
				a.chrome = "24";
			
			if ( a.chrome.compare("23") == 0 )
				a.chrome = "chrX";
			else if ( a.chrome.compare("24") == 0 )
				a.chrome = "chrY";
	        else
				a.chrome = string("chr") + a.chrome;
			
			if(data.find(a.target) == data.end()) {
				map<string, list<segment> > l;
				data[ a.target ] = l;
			} 
			if ( data[ a.target ].find( a.chrome ) == data[ a.target ].end() ) {
				list<segment> l;
				data[ a.target ][ a.chrome ] = l;
			}
			data[ a.target ][a.chrome].push_back( a );
			targetSet.insert( a.target );
		}
	} while ( cin );	
	
	//create a break map
	breakmap breaks;
	for( segmap::iterator t = data.begin(); t!=data.end(); ++t) {
		for ( chromemap::iterator c = t->second.begin(); c != t->second.end(); ++c) {
			for ( list<segment>::iterator s = c->second.begin(); s != c->second.end(); s++ ) {
				if ( breaks.find( s->chrome ) == breaks.end() ) {
					set<int> ns;
					ns.insert( s->start );
					ns.insert( s->end );
					breaks[ s->chrome ] = ns;
				} else {
					breaks[ s->chrome ].insert( s->start );
					breaks[ s->chrome ].insert( s->end );
				}
			}
		}	
	}
	
	//print out the matrix column names
	cout << "probe";
	for ( set<string>::iterator ts = targetSet.begin(); ts != targetSet.end(); ++ts ) {
		cout << "\t";
		cout << (*ts);
	}
	cout << "\n";
	int targetCount = targetSet.size();
	//check breakmap for each chrome
	for ( breakmap::iterator b = breaks.begin(); b != breaks.end(); ++b ) {
		cerr << b->first << "\n";
		vector<int> starts;
		vector<int> ends;
		vector<string> probeName;
		
		breakmap::iterator cset = breaks.find( b->first );
		int blockCount = cset->second.size() - 1;
		starts.resize( blockCount );
		ends.resize( blockCount );
		probeName.resize( blockCount );
		int i = 0;
		for ( set<int>::iterator v = cset->second.begin(); v != cset->second.end(); ++v ) {
			if ( i > 0 ) {
				ends[i-1] = (*v);
			}
			if ( i < blockCount ) {
				starts[ i ] = (*v);
			}
			i++;
		}
		
		//create names for the breakpoints
		genemap gm;
		gm.resize( blockCount );
		for ( int i = 0; i < blockCount; i++ ) {
			stringstream name;
			name << b->first;
			name << "_";
			name << starts[i];
			name << "_";
			name << ends[i];
			probeName[i] = name.str();
			gm[ i ].resize( targetCount ); 
		}
		//scan the targets		
		int curTarget = 0;
		for ( set<string>::iterator t = targetSet.begin(); t != targetSet.end(); ++t ) {
			//find the target values for the current chrome
			segmap::iterator s = data.find( (*t) );
			chromemap::iterator ss = s->second.find( b->first );
			if ( ss != s->second.end() ) {
				vector<float> vals;
				vals.resize( blockCount );
				for ( int i = 0; i < blockCount; i++ ) {
					vals[ i ] = MISSING_VAL;
				}
				//if the segment overlaps the break, assign the value
				for ( list<segment>::iterator cs = ss->second.begin(); cs != ss->second.end(); cs++ ) {
					for ( int i = 0; i < blockCount; i++ ) {
						if ( cs->end > starts[i] && cs->start < ends[i] ) {
							vals[i] = cs->value;
						}
					}
				}
				//assign the values to the named map
				for ( int i = 0; i < blockCount; i++ ) {
					gm[ i ][ curTarget ] = vals[ i ];
				}				
				//cout << s->first << "\t" << ss->first << "\n";
			}
			curTarget++;
		}
		//print out this chrome's segments
		//and the values for each target
		for ( int i =0; i < blockCount; i++ ) {
			cout << probeName[i];
			for ( int j = 0; j < targetCount; j++ ) {
				float val = gm[i][j];
				cout << "\t";
				if ( val == MISSING_VAL )
					cout << "null";
				else
					cout << val;
			}
			cout << "\n";
		}
	}
	
	return 0;
}

#ifdef __cplusplus
}
#endif
