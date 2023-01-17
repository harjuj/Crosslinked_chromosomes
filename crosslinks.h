#ifndef CROSSLINKS_H
#define CROSSLINKS_H

//Given two REDUCED indices, remove the crosslink between them
void detach(int site, int other, int thread_num){
    auto itj = std::find(crosslinked_with[thread_num][other].begin(), crosslinked_with[thread_num][other].end(), site);
    if (itj==crosslinked_with[thread_num][other].end()){
        std::cerr<<"Error: didn't find other site's crosslink: ("<<site<<")..."<<std::endl;
        exit(EXIT_FAILURE);
    }
    crosslinked_with[thread_num][other].erase(itj);
    if (crosslinked_with[thread_num][other].size()==0){
        crosslinked_with[thread_num].erase(other);        
        is_crosslinked[thread_num][other]=false;        
    }
    
    auto iti = std::find(crosslinked_with[thread_num][site].begin(), crosslinked_with[thread_num][site].end(), other);
    if (iti==crosslinked_with[thread_num][site].end()){
        std::cerr<<"Error: didn't find site's crosslink: ("<<other<<")..."<<std::endl;
        exit(EXIT_FAILURE);
    }
    crosslinked_with[thread_num][site].erase(iti);
    if (crosslinked_with[thread_num][site].size()==0){
        crosslinked_with[thread_num].erase(site);
        is_crosslinked[thread_num][site]=false;        
    }
    current_crosslinks[thread_num]--;
}

//Given monomer index, try detaching each of its crosslinks
void try_detach(int mon, int thread_num){
    //make a list of monomers to be detached
    if (mon%reduction_factor==0 && is_crosslinked[thread_num][mon/reduction_factor]){
        std::vector<int> to_detach;
        int site{mon/reduction_factor};

        for (int other : crosslinked_with[thread_num][site]){
            if (generators[thread_num].disReal()<std::min(1.,p_det*exp(binding_affinities[site]+binding_affinities[other]))){
                to_detach.push_back(other);
            }
        }    
        for (int other : to_detach){
            detach(site, other, thread_num);
        }    
    }
}

//Given REDUCED monomer indices, attach crosslink between them
void attach(int site, int other, int thread_num){
    if (is_crosslinked[thread_num][site]){
        crosslinked_with[thread_num][site].push_back(other);    
    }
    else{
        crosslinked_with[thread_num][site]={other};
        is_crosslinked[thread_num][site]=true;
    }

    if (is_crosslinked[thread_num][other]){
        crosslinked_with[thread_num][other].push_back(site);    
    }
    else{
        crosslinked_with[thread_num][other]={site};
        is_crosslinked[thread_num][other]=true;
    }
    current_crosslinks[thread_num]++;
}

//Given monomer index, for each contact try attaching a crosslink
void try_attach(int mon, int thread_num){
    if (mon%reduction_factor==0){
        std::vector<int> x{polymer[thread_num][mon][0], polymer[thread_num][mon][1], polymer[thread_num][mon][2]};
        
        int site=mon/reduction_factor;

        for (int other : locations[thread_num][x]){
            if (other!=site && generators[thread_num].disReal()<std::min(1.,p_cl*exp(-binding_affinities[site]-binding_affinities[other]))){
                attach(site, other, thread_num);
            }
        }    
    }    
}

#endif
