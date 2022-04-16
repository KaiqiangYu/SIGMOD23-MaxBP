#pragma once
#include<iostream>
#include<vector>
#include<time.h>
#include"RandList.h"
#include"subgraph.h"
#include<queue>
using namespace std;

class BKMB
{
private:
    /* data */
    int **Graph;
    int *degree;
    int graph_size;
    int bi_index;
    int k;
    int left, right;
    int temp_i,temp_j,temp_k, temp_value, temp_node, temp_node2, temp_size1, temp_size2;
    MBitSet *bit_G;
    RandList P_L, P_R, Cand_L, Cand_R, Exc_L, Exc_R;
    int *degInP;
    int *degInG;
    int *G_index; int count;
    vector<int> temp_vector;
    int topK;
    bool CanBranch();
    void BKmb_Rec();
    void Record(vector<int>& left, vector<int>& right);

public:
    BKMB(int **Graph, int *degree, int graph_size, int bi_index, int k,int left, int right, int topK, priority_queue<subgraph,vector<subgraph>,greater<subgraph>> top_k);
    ~BKMB();
    void BKmb();
    int total_num;
    int tt;
    int edges,v_l,v_r;
    int myset;
    //#ifdef OUTPUT
    int *re_temp_index;
    // subgraph *maxbp;
    priority_queue<subgraph,vector<subgraph>,greater<subgraph>> top_k;
    //#endif
};

BKMB::BKMB(int **Graph, int *degree, int graph_size, int bi_index, int k, int left, int right, int topK, priority_queue<subgraph,vector<subgraph>,greater<subgraph>> top_k)
{
    this->Graph=Graph;
    this->degree=degree;
    this->graph_size=graph_size;
    this->bi_index=bi_index;
    this->k=k;
    this->left=left;
    this->right=right;
    total_num=0;
    edges=0;
    myset=graph_size;

    P_L.init(graph_size);
    P_R.init(graph_size);
    Cand_L.init(graph_size);
    Cand_R.init(graph_size);
    Exc_L.init(graph_size);
    Exc_R.init(graph_size);
    temp_vector.reserve(graph_size);
    tt=0;
    this->top_k=top_k;
    this->topK=topK;
    if(top_k.empty()) edges=0;
    else{
        subgraph temp=top_k.top();
        edges=temp.edges;
    }

    degInG=new int[graph_size];
    degInP=new int[graph_size];
    G_index=new int[graph_size];
    count=0;

    for(temp_i=graph_size-1;temp_i>=0;--temp_i){
        degInP[temp_i]=0;
        degInG[temp_i]=0;
        G_index[temp_i]=0;
    }
    

    bit_G=new MBitSet[graph_size];
    double ttm=0;
    for(int i=0;i<graph_size;++i){
        bit_G[i].allocacte(graph_size+1);
        bit_G[i].reinit(graph_size+1);
        for(int j=0;j<degree[i];++j){
            bit_G[i].set(Graph[i][j]);
        }
        ttm+=sizeof(bit_G[0]);
    }
}

BKMB::~BKMB()
{
    delete [] degInG;
    delete [] degInP;
    delete [] G_index;
    //delete []bit_G;
    //delete &temp_vector;
    //delete &P_L, &P_R, &Cand_L, &Cand_R, &Exc_L, &Exc_R;
}

void BKMB::Record(vector<int>& left, vector<int>& right){
    for(int i=P_L.vnum-1;i>=0;--i)
        left.push_back(re_temp_index[P_L.vlist[i]]);
    for(int i=Cand_L.vnum-1;i>=0;--i)
        left.push_back(re_temp_index[Cand_L.vlist[i]]);
    for(int i=P_R.vnum-1;i>=0;--i)
        right.push_back(re_temp_index[P_R.vlist[i]]);
    for(int i=Cand_R.vnum-1;i>=0;--i)
        right.push_back(re_temp_index[Cand_R.vlist[i]]);
}

void BKMB::BKmb(){
    bool flag;
    for(int i=0;i<bi_index;++i){
        P_L.add(i);
        for(temp_i=degree[i]-1;temp_i>=0;--temp_i){
            temp_node=Graph[i][temp_i];
            count=0;
            for(temp_j=degree[temp_node]-1;temp_j>=0;--temp_j){
                temp_node2=Graph[temp_node][temp_j];
                if(temp_node2>=i) count++;
            }
            if(count<left-k) continue;
            degInG[temp_node]++;
            degInP[temp_node]++;
            Cand_R.add(temp_node);
            for(temp_j=degree[temp_node]-1;temp_j>=0;--temp_j){
                temp_node2=Graph[temp_node][temp_j];
                degInG[temp_node2]++;
            }
        }

        for(temp_i=degree[i]-1;temp_i>=0;--temp_i){
            temp_node=Graph[i][temp_i];
            for(temp_j=degree[temp_node]-1;temp_j>=0;--temp_j){
                temp_node2=Graph[temp_node][temp_j];
                if(temp_node2>i&&degInG[temp_node2]>=right-2*k&&Cand_L.vpos[temp_node2]==graph_size){
                    Cand_L.add(temp_node2);
                    for(temp_k=degree[temp_node2]-1;temp_k>=0;--temp_k){
                        degInG[Graph[temp_node2][temp_k]]++;
                    }
                }
                if(temp_node2<i&&degInG[temp_node2]>=right-2*k&&Exc_L.vpos[temp_node2]==graph_size){
                    Exc_L.add(temp_node2);
                }
            }
        }

        for(temp_i=bi_index;temp_i<graph_size;++temp_i){
            if(degInG[temp_i]>=left-k&&Cand_R.vpos[temp_i]==graph_size){
                Cand_R.add(temp_i);
                for(temp_j=degree[temp_i]-1;temp_j>=0;--temp_j){
                    degInG[Graph[temp_i][temp_j]]++;
                }
            }
            if(degInG[temp_i]<left-k&&Cand_R.vpos[temp_i]<graph_size){
                Cand_R.remove(temp_i);
                for(temp_j=degree[temp_i]-1;temp_j>=0;--temp_j){
                    temp_node=Graph[temp_i][temp_j];
                    degInG[temp_node]--;
                    if(degInG[temp_node]<right-k&&Cand_L.vpos[temp_node]<graph_size){
                        Cand_L.remove(temp_node);
                        for(temp_k=degree[temp_node]-1;temp_k>=0;--temp_k){
                            --degInG[Graph[temp_node][temp_k]];
                        }
                    }
                }
            }
        }

        for(temp_i=bi_index;temp_i<graph_size;++temp_i){
            if(degInG[temp_i]>=left-k&&Cand_R.vpos[temp_i]==graph_size){
                Cand_R.add(temp_i);
                for(temp_j=degree[temp_i]-1;temp_j>=0;--temp_j){
                    degInG[Graph[temp_i][temp_j]]++;
                }
            }
            if(degInG[temp_i]<left-k&&Cand_R.vpos[temp_i]<graph_size){
                Cand_R.remove(temp_i);
                for(temp_j=degree[temp_i]-1;temp_j>=0;--temp_j){
                    temp_node=Graph[temp_i][temp_j];
                    degInG[temp_node]--;
                    if(degInG[temp_node]<right-k&&Cand_L.vpos[temp_node]<graph_size){
                        Cand_L.remove(temp_node);
                        for(temp_k=degree[temp_node]-1;temp_k>=0;--temp_k){
                            --degInG[Graph[temp_node][temp_k]];
                        }
                    }
                }
            }
        }

        for(temp_i=bi_index;temp_i<graph_size;++temp_i){
            if(degInG[temp_i]>=left-k&&Cand_R.vpos[temp_i]==graph_size){
                Cand_R.add(temp_i);
                for(temp_j=degree[temp_i]-1;temp_j>=0;--temp_j){
                    degInG[Graph[temp_i][temp_j]]++;
                }
            }
            if(degInG[temp_i]<left-k&&Cand_R.vpos[temp_i]<graph_size){
                Cand_R.remove(temp_i);
                for(temp_j=degree[temp_i]-1;temp_j>=0;--temp_j){
                    temp_node=Graph[temp_i][temp_j];
                    degInG[temp_node]--;
                    if(degInG[temp_node]<right-k&&Cand_L.vpos[temp_node]<graph_size){
                        Cand_L.remove(temp_node);
                        for(temp_k=degree[temp_node]-1;temp_k>=0;--temp_k){
                            --degInG[Graph[temp_node][temp_k]];
                        }
                    }
                }
            }
        }

        flag=true;
        for(temp_i=Exc_L.vnum-1;temp_i>=0;--temp_i){
            temp_node=Exc_L.vlist[temp_i];
            if(degInG[temp_node]==Cand_R.vnum){
                flag=false;
                break;
            }
        }
        
        if(flag)
        BKmb_Rec();

        for(temp_i=Cand_R.vnum-1;temp_i>=0;--temp_i){
            temp_node=Cand_R.vlist[temp_i];
            for(temp_j=degree[temp_node]-1;temp_j>=0;--temp_j){
                --degInG[Graph[temp_node][temp_j]];
            }
        }
        for(temp_i=Cand_L.vnum-1;temp_i>=0;--temp_i){
            temp_node=Cand_L.vlist[temp_i];
            for(temp_j=degree[temp_node]-1;temp_j>=0;--temp_j){
                --degInG[Graph[temp_node][temp_j]];
            }
        }
        
        for(temp_i=P_L.vnum-1;temp_i>=0;--temp_i){
            temp_node=P_L.vlist[temp_i];
            for(temp_j=degree[temp_node]-1;temp_j>=0;--temp_j){
                --degInG[Graph[temp_node][temp_j]];
                --degInP[Graph[temp_node][temp_j]];
            }
        }
        
        Cand_R.clear();
        Cand_L.clear();
        P_L.clear();
        Exc_L.clear();      
    }
}

void BKMB::BKmb_Rec(){
    if(P_L.vnum>=myset){
        return;
    }
    
    if(P_L.vnum+Cand_L.vnum<left||P_R.vnum+Cand_R.vnum<right) return;
    
    if((P_L.vnum+Cand_L.vnum)*(P_R.vnum+Cand_R.vnum)<=edges&&(P_L.vnum+Cand_L.vnum)*(P_R.vnum+Cand_R.vnum)>0) return;
    
    temp_j=0;
    for(temp_i=P_L.vnum-1;temp_i>=0;--temp_i){
        temp_node=P_L.vlist[temp_i];
        temp_j+=degInG[temp_node];
    }
    for(temp_i=Cand_L.vnum-1;temp_i>=0;--temp_i){
        temp_node=Cand_L.vlist[temp_i];
        temp_j+=degInG[temp_node];
    }
    if(temp_j<=edges){
        return;
    }
    
    if(Cand_R.empty()&&Cand_L.empty()){
        if(Exc_R.empty()&&Exc_L.empty()){
            total_num++;
            temp_j=0;
            for(temp_i=P_L.vnum-1;temp_i>=0;--temp_i){
                temp_node=P_L.vlist[temp_i];
                temp_j+=degInP[temp_node];
            }
            if(edges<temp_j){
                // edges=temp_j;
                // v_l=P_L.vnum+Cand_L.vnum;
                // v_r=P_R.vnum+Cand_R.vnum;

                //#ifdef OUTPUT
                // maxbp->edges=temp_j;
                // maxbp->left_size=v_l;
                // maxbp->right_size=v_r;
                // maxbp->left.clear();
                // maxbp->right.clear();
                // for(int ptr_i=P_L.vnum-1;ptr_i>=0;--ptr_i)  maxbp->left.push_back(re_temp_index[P_L.vlist[ptr_i]]);
                // for(int ptr_i=Cand_L.vnum-1;ptr_i>=0;--ptr_i)  maxbp->left.push_back(re_temp_index[Cand_L.vlist[ptr_i]]);
                // for(int ptr_i=P_R.vnum-1;ptr_i>=0;--ptr_i)  maxbp->right.push_back(re_temp_index[P_R.vlist[ptr_i]]);
                // for(int ptr_i=Cand_R.vnum-1;ptr_i>=0;--ptr_i)  maxbp->right.push_back(re_temp_index[Cand_R.vlist[ptr_i]]);
                //#endif
                
                if(top_k.size()<topK){
                    subgraph temp_sub(temp_j,P_L.vnum+Cand_L.vnum,P_R.vnum+Cand_R.vnum);
                    Record(temp_sub.left,temp_sub.right);
                    top_k.push(temp_sub);
                    //edges=edges<temp_j&&edges>0?edges:temp_j;
                }else if(edges<temp_j){
                    subgraph temp_sub(temp_j,P_L.vnum+Cand_L.vnum,P_R.vnum+Cand_R.vnum);
                    Record(temp_sub.left,temp_sub.right);
                    top_k.push(temp_sub);
                    top_k.pop();
                    edges=top_k.top().edges;
                }
            }

        } 
        return;
    }

    int povit=-1;
    temp_value=0;
    temp_node2=graph_size;
    temp_j=P_L.vnum+Cand_L.vnum;
    temp_k=P_R.vnum+Cand_R.vnum;

    for(temp_i=P_L.vnum-1;temp_i>=0;--temp_i){
        temp_node=P_L.vlist[temp_i];
        if(temp_k-degInG[temp_node]>temp_value){
            temp_value=temp_k-degInG[temp_node];
            povit=temp_node;
        }
        if(temp_node2>degInG[temp_node]){
            temp_node2=degInG[temp_node];
        }
    }
    if(temp_node2+k<right){
        return;
    }
    temp_size1=graph_size;
    for(temp_i=P_R.vnum-1;temp_i>=0;--temp_i){
        temp_node=P_R.vlist[temp_i];
        if(temp_j-degInG[temp_node]>temp_value){
            temp_value=temp_j-degInG[temp_node];
            povit=temp_node;
        }
        if(temp_size1>degInG[temp_node]){
            temp_size1=degInG[temp_node];
        }
    }
    if(temp_size1+k<left||(temp_node2+k)*(temp_size1+k)<=edges+2*k){
        return;
    }
    
    if(povit>=0&&temp_value>k){
        RandList *Cand1, *Cand2, *P1, *P2, *Exc1, *Exc2;
        int temp_stand1, temp_stand2;
        if(povit<bi_index){
            Cand1=&Cand_L;
            Cand2=&Cand_R;
            P1=&P_L;
            P2=&P_R;
            Exc1=&Exc_L;
            Exc2=&Exc_R;
            temp_stand1=left;
            temp_stand2=right;
        }else{
            Cand1=&Cand_R;
            Cand2=&Cand_L;
            P1=&P_R;
            P2=&P_L;
            Exc1=&Exc_R;
            Exc2=&Exc_L;
            temp_stand1=right;
            temp_stand2=left;
        }

        vector<int> doing;
        doing.reserve(Cand2->vnum);
        for(temp_i=Cand2->vnum-1;temp_i>=0;--temp_i){
            temp_node=Cand2->vlist[temp_i];
            if(!bit_G[povit].test(temp_node)){
                doing.push_back(temp_node);
            }
        }
        int p=k-(P2->vnum-degInP[povit]), idx=1;
        int rec_a, rec_b;
        rec_a=doing[0];
        Cand2->remove(rec_a);
        Exc2->add(rec_a);
        
        for(temp_i=degree[rec_a]-1;temp_i>=0;--temp_i){
            --degInG[Graph[rec_a][temp_i]];
        }
        if(CanBranch()&&P2->vnum+Cand2->vnum>=temp_stand2)
            BKmb_Rec();

        Exc2->remove(rec_a);
        vector<int> remove_C1, remove_E1, remove_C2, remove_E2;
        remove_C1.reserve(Cand1->vnum);
        remove_C2.reserve(Cand2->vnum);
        remove_E1.reserve(Exc1->vnum);
        remove_E2.reserve(Exc2->vnum);
      
        while(idx<p+1){
            
            P2->add(rec_a);
            for(temp_i=degree[rec_a]-1;temp_i>=0;--temp_i){
                temp_node=Graph[rec_a][temp_i];
                degInG[temp_node]++;
                degInP[temp_node]++;
            }

            if(degInP[rec_a]==P1->vnum-k){
                for(temp_i=Cand1->vnum-1;temp_i>=0;--temp_i){
                    temp_node=Cand1->vlist[temp_i];
                    if(!bit_G[rec_a].test(temp_node)||degInG[temp_node]<temp_stand2-k){
                        Cand1->remove(temp_node);
                        remove_C1.push_back(temp_node);
                        for(temp_j=degree[temp_node]-1;temp_j>=0;--temp_j){
                            temp_node2=Graph[temp_node][temp_j];
                            --degInG[temp_node2];
                        }
                    }
                    
                }

                for(temp_i=Exc1->vnum-1;temp_i>=0;--temp_i){
                    temp_node=Exc1->vlist[temp_i];
                    if(!bit_G[rec_a].test(temp_node)||degInG[temp_node]<temp_stand2-k){
                        Exc1->remove(temp_node);
                        remove_E1.push_back(temp_node);
                    }
                }
            }else{
                temp_size1=P2->vnum-k;
                for(temp_i=Cand1->vnum-1;temp_i>=0;--temp_i){
                    temp_node=Cand1->vlist[temp_i];
                    if(degInP[temp_node]<temp_size1||degInG[temp_node]<temp_stand2-k){
                        Cand1->remove(temp_node);
                        remove_C1.push_back(temp_node);
                        for(temp_j=degree[temp_node]-1;temp_j>=0;--temp_j){
                            temp_node2=Graph[temp_node][temp_j];
                            --degInG[temp_node2];
                        }
                    }
                }

                for(temp_i=Exc1->vnum-1;temp_i>=0;--temp_i){
                    temp_node=Exc1->vlist[temp_i];
                    if(degInP[temp_node]<temp_size1||degInG[temp_node]<temp_stand2-k){
                        Exc1->remove(temp_node);
                        remove_E1.push_back(temp_node);
                    }
                }
            }
            
            if(Cand1->vnum+P1->vnum<temp_stand1){
                break;
            }

            temp_vector.clear();
            temp_size1=P2->vnum;
            temp_size2=P1->vnum;
            for(temp_i=P1->vnum-1;temp_i>=0;--temp_i){
                temp_node=P1->vlist[temp_i];
                if(degInP[temp_node]==temp_size1-k&&!bit_G[rec_a].test(temp_node)){
                    temp_vector.push_back(temp_node);
                }
            }
            for(temp_i=degree[rec_a]-1;temp_i>=0;--temp_i){
                temp_node=Graph[rec_a][temp_i];
                if(Cand1->contains(temp_node)||P1->contains(temp_node))
                    G_index[temp_node]=1;
            }
                    
            
            if(!temp_vector.empty()){
                for(temp_i=Cand2->vnum-1;temp_i>=0;--temp_i){
                    temp_node=Cand2->vlist[temp_i];
                    if(degInP[temp_node]==temp_size2) continue;
                    if(degInP[temp_node]<temp_size2-k||degInG[temp_node]<temp_stand1-k){
                        Cand2->remove(temp_node);
                        remove_C2.push_back(temp_node);
                        for(temp_j=degree[temp_node]-1;temp_j>=0;--temp_j){
                            temp_node2=Graph[temp_node][temp_j];
                            --degInG[temp_node2];
                        }
                    }else{
                        for(int i:temp_vector){
                            if(!bit_G[temp_node].test(i)){
                                Cand2->remove(temp_node);
                                remove_C2.push_back(temp_node);
                                for(temp_j=degree[temp_node]-1;temp_j>=0;--temp_j){
                                    temp_node2=Graph[temp_node][temp_j];
                                    --degInG[temp_node2];
                                }
                                break;
                            }
                        }

                        if(Cand2->contains(temp_node)){
                            count=0;
                            for(temp_j=degree[temp_node]-1;temp_j>=0;--temp_j){
                                if(G_index[Graph[temp_node][temp_j]]){
                                    count++;
                                }
                            }
                            if(count<temp_stand1-2*k){
                                Cand2->remove(temp_node);
                                remove_C2.push_back(temp_node);
                                for(temp_j=degree[temp_node]-1;temp_j>=0;--temp_j){
                                    temp_node2=Graph[temp_node][temp_j];
                                    --degInG[temp_node2];
                                }
                            }
                        }   
                    }
                }

                for(temp_i=Exc2->vnum-1;temp_i>=0;--temp_i){
                    temp_node=Exc2->vlist[temp_i];
                    if(degInP[temp_node]==temp_size2) continue;
                    if(degInP[temp_node]<temp_size2-k||degInG[temp_node]<temp_stand1-k){
                        Exc2->remove(temp_node);
                        remove_E2.push_back(temp_node);
                    }else{
                        for(int i:temp_vector){
                            if(!bit_G[temp_node].test(i)){
                                Exc2->remove(temp_node);
                                remove_E2.push_back(temp_node);
                                break;
                            }
                        }
                    }
                }
            }

            for(temp_i=degree[rec_a]-1;temp_i>=0;--temp_i){
                temp_node=Graph[rec_a][temp_i];
                G_index[temp_node]=0;
            }
            

            if(Cand2->vnum+P2->vnum<temp_stand2){
                break;
            }

            for(;idx<doing.size()&&!Cand2->contains(doing[idx]);++p,++idx);
            if(idx==doing.size()) break;
            rec_b=doing[idx];

            Cand2->remove(rec_b);
            Exc2->add(rec_b);
            for(temp_i=degree[rec_b]-1;temp_i>=0;--temp_i){
                temp_node=Graph[rec_b][temp_i];
                --degInG[temp_node];
            }
            if(CanBranch())
                BKmb_Rec();

            Exc2->remove(rec_b);
            rec_a=rec_b;
            ++idx;
        }

        if(Cand1->vnum+P1->vnum>=temp_stand1&&Cand2->vnum+P2->vnum>=temp_stand2){
            
            count=0;
            for(;idx<doing.size();++idx){
                temp_node=doing[idx];
                if(Cand2->contains(temp_node)){
                    count++;
                }
            }
            
            if(Cand2->vnum+P2->vnum-count>=temp_stand2){
                for(;idx<doing.size();++idx){
                temp_node=doing[idx];
                if(Cand2->contains(temp_node)){
                    Cand2->remove(temp_node);
                    for(temp_i=degree[temp_node]-1;temp_i>=0;--temp_i){
                        --degInG[Graph[temp_node][temp_i]];
                    }
                }
                }
                if(CanBranch())
                    BKmb_Rec();
            }
            
            
            
        }
        


        
        for(int i:doing){
            if(P2->vpos[i]<graph_size){
                P2->remove(i);
                Cand2->add(i);
                for(temp_i=degree[i]-1;temp_i>=0;--temp_i){
                    temp_node=Graph[i][temp_i];
                    --degInP[temp_node];
                }
            }
            if(Cand2->vpos[i]>=graph_size){
                Cand2->add(i);
                for(temp_i=degree[i]-1;temp_i>=0;--temp_i){
                    temp_node=Graph[i][temp_i];
                    degInG[temp_node]++;
                }

            }
        }

        for(int i:remove_C2){
            if(Cand2->vpos[i]==graph_size){
                Cand2->add(i);
                for(temp_i=degree[i]-1;temp_i>=0;--temp_i){
                    temp_node=Graph[i][temp_i];
                    degInG[temp_node]++;
                }
            }
        }

        for(int i:remove_C1){
            Cand1->add(i);
            for(temp_i=degree[i]-1;temp_i>=0;--temp_i){
                temp_node=Graph[i][temp_i];
                degInG[temp_node]++;
            }
        }

        for(int i:remove_E1){
            Exc1->add(i);
        }
        for(int i:remove_E2){
            Exc2->add(i);
        }

        

    }else{
        for(temp_i=Cand_L.vnum-1;temp_i>=0;--temp_i){
            temp_node=Cand_L.vlist[temp_i];
            if(temp_k-degInG[temp_node]>temp_value){
                temp_value=temp_k-degInG[temp_node];
                povit=temp_node;
            }
        }

        for(temp_i=Cand_R.vnum-1;temp_i>=0;--temp_i){
            temp_node=Cand_R.vlist[temp_i];
            if(temp_j-degInG[temp_node]>temp_value){
                temp_value=temp_j-degInG[temp_node];
                povit=temp_node;
            }
        }

        if(temp_value<=k){
            if(Exc_L.empty()&&Exc_R.empty()){
                ++total_num;

                temp_value=0;
                for(temp_i=P_L.vnum-1;temp_i>=0;--temp_i){
                    temp_node=P_L.vlist[temp_i];
                    temp_value+=degInG[temp_node];
                }
                for(temp_i=Cand_L.vnum-1;temp_i>=0;--temp_i){
                    temp_node=Cand_L.vlist[temp_i];
                    temp_value+=degInG[temp_node];
                }
                if(edges<temp_value){
                    // edges=temp_value;
                    // v_l=P_L.vnum+Cand_L.vnum;
                    // v_r=P_R.vnum+Cand_R.vnum;


                    //#ifdef OUTPUT
                    // maxbp->edges=temp_value;
                    // maxbp->left_size=v_l;
                    // maxbp->right_size=v_r;
                    // maxbp->left.clear();
                    // maxbp->right.clear();
                    // for(int ptr_i=P_L.vnum-1;ptr_i>=0;--ptr_i)  maxbp->left.push_back(re_temp_index[P_L.vlist[ptr_i]]);
                    // for(int ptr_i=Cand_L.vnum-1;ptr_i>=0;--ptr_i)  maxbp->left.push_back(re_temp_index[Cand_L.vlist[ptr_i]]);
                    // for(int ptr_i=P_R.vnum-1;ptr_i>=0;--ptr_i)  maxbp->right.push_back(re_temp_index[P_R.vlist[ptr_i]]);
                    // for(int ptr_i=Cand_R.vnum-1;ptr_i>=0;--ptr_i)  maxbp->right.push_back(re_temp_index[Cand_R.vlist[ptr_i]]);
                    //#endif
                    
                    if(top_k.size()<topK){
                        subgraph temp_sub(temp_value,P_L.vnum+Cand_L.vnum,P_R.vnum+Cand_R.vnum);
                        Record(temp_sub.left,temp_sub.right);
                        top_k.push(temp_sub);
                    //edges=edges<temp_j&&edges>0?edges:temp_j;
                    }else if(edges<temp_value){
                        subgraph temp_sub(temp_value,P_L.vnum+Cand_L.vnum,P_R.vnum+Cand_R.vnum);
                        Record(temp_sub.left,temp_sub.right);
                        top_k.push(temp_sub);
                        top_k.pop();
                        edges=top_k.top().edges;
                    }
                }
                return;
            }

            for(temp_i=Exc_R.vnum-1;temp_i>=0;--temp_i){
                temp_node=Exc_R.vlist[temp_i];
                if(degInG[temp_node]==temp_j) return;
            }

            for(temp_i=Exc_L.vnum-1;temp_i>=0;--temp_i){
                temp_node=Exc_L.vlist[temp_i];
                if(degInG[temp_node]==temp_k) return;
            }

            temp_vector.clear();
            for(temp_i=P_L.vnum-1;temp_i>=0;--temp_i){
                temp_node=P_L.vlist[temp_i];
                if(degInG[temp_node]==temp_k-k){
                    temp_vector.push_back(temp_node);
                }
            }
            for(temp_i=Cand_L.vnum-1;temp_i>=0;--temp_i){
                temp_node=Cand_L.vlist[temp_i];
                if(degInG[temp_node]==temp_k-k){
                    temp_vector.push_back(temp_node);
                }
            }
            for(temp_i=Exc_R.vnum-1;temp_i>=0;--temp_i){
                temp_node=Exc_R.vlist[temp_i];
                if(degInG[temp_node]>=temp_j-k){
                    temp_value=1;
                    for(int i:temp_vector){
                        if(!bit_G[temp_node].test(i)){
                            temp_value=0;
                            break;
                        }
                    }
                    if(temp_value) return;
                }
            }


            temp_vector.clear();
            for(temp_i=P_R.vnum-1;temp_i>=0;--temp_i){
                temp_node=P_R.vlist[temp_i];
                if(degInG[temp_node]==temp_j-k){
                    temp_vector.push_back(temp_node);
                }
            }
            for(temp_i=Cand_R.vnum-1;temp_i>=0;--temp_i){
                temp_node=Cand_R.vlist[temp_i];
                if(degInG[temp_node]==temp_j-k){
                    temp_vector.push_back(temp_node);
                }
            }
            for(temp_i=Exc_L.vnum-1;temp_i>=0;--temp_i){
                temp_node=Exc_L.vlist[temp_i];
                if(degInG[temp_node]>=temp_k-k){
                    temp_value=1;
                    for(int i:temp_vector){
                        if(!bit_G[temp_node].test(i)){
                            temp_value=0;
                            break;
                        }
                    }
                    if(temp_value) return;
                }
            }


            temp_value=0;
            for(temp_i=P_L.vnum-1;temp_i>=0;--temp_i){
                temp_node=P_L.vlist[temp_i];
                temp_value+=degInG[temp_node];
            }
            for(temp_i=Cand_L.vnum-1;temp_i>=0;--temp_i){
                temp_node=Cand_L.vlist[temp_i];
                temp_value+=degInG[temp_node];
            }
            if(edges<temp_value){
                // edges=temp_value;
                // v_l=P_L.vnum+Cand_L.vnum;
                // v_r=P_R.vnum+Cand_R.vnum;

                //#ifdef OUTPUT
                // maxbp->edges=temp_value;
                // maxbp->left_size=v_l;
                // maxbp->right_size=v_r;
                // maxbp->left.clear();
                // maxbp->right.clear();
                // for(int ptr_i=P_L.vnum-1;ptr_i>=0;--ptr_i)  maxbp->left.push_back(re_temp_index[P_L.vlist[ptr_i]]);
                // for(int ptr_i=Cand_L.vnum-1;ptr_i>=0;--ptr_i)  maxbp->left.push_back(re_temp_index[Cand_L.vlist[ptr_i]]);
                // for(int ptr_i=P_R.vnum-1;ptr_i>=0;--ptr_i)  maxbp->right.push_back(re_temp_index[P_R.vlist[ptr_i]]);
                // for(int ptr_i=Cand_R.vnum-1;ptr_i>=0;--ptr_i)  maxbp->right.push_back(re_temp_index[Cand_R.vlist[ptr_i]]);
                //#endif
                
                if(top_k.size()<topK){
                    subgraph temp_sub(temp_value,P_L.vnum+Cand_L.vnum,P_R.vnum+Cand_R.vnum);
                    Record(temp_sub.left,temp_sub.right);
                    top_k.push(temp_sub);
                    //edges=edges<temp_value&&edges>0?edges:temp_value;
                }else if(edges<temp_value){
                    subgraph temp_sub(temp_value,P_L.vnum+Cand_L.vnum,P_R.vnum+Cand_R.vnum);
                    Record(temp_sub.left,temp_sub.right);
                    top_k.push(temp_sub);
                    top_k.pop();
                    edges=top_k.top().edges;
                }
            }

            ++total_num;
            return;
        }

        
        RandList *Cand1, *Cand2, *P1, *P2, *Exc1, *Exc2;
        int temp_stand1, temp_stand2;
        if(povit<bi_index){
            Cand1=&Cand_L;
            Cand2=&Cand_R;
            P1=&P_L;
            P2=&P_R;
            Exc1=&Exc_L;
            Exc2=&Exc_R;
            temp_stand1=left;
            temp_stand2=right;
        }else{
            Cand1=&Cand_R;
            Cand2=&Cand_L;
            P1=&P_R;
            P2=&P_L;
            Exc1=&Exc_R;
            Exc2=&Exc_L;
            temp_stand1=right;
            temp_stand2=left;
        }

        Cand1->remove(povit);
        Exc1->add(povit);
        for(temp_i=degree[povit]-1;temp_i>=0;--temp_i){
            temp_node=Graph[povit][temp_i];
            --degInG[temp_node];
        }
        
        if(CanBranch())
            BKmb_Rec();

        Exc1->remove(povit);
        P1->add(povit);
        for(temp_i=degree[povit]-1;temp_i>=0;--temp_i){
            temp_node=Graph[povit][temp_i];
            degInG[temp_node]++;
            degInP[temp_node]++;
        }

        
        vector<int> remove_C1, remove_E1, remove_C2, remove_E2;
        remove_C1.reserve(Cand1->vnum);
        remove_C2.reserve(Cand2->vnum);
        remove_E1.reserve(Exc1->vnum);
        remove_E2.reserve(Exc2->vnum);

        
        
        if(degInP[povit]==P2->vnum-k){
            for(temp_i=Cand2->vnum-1;temp_i>=0;--temp_i){
                temp_node=Cand2->vlist[temp_i];
                if(!bit_G[povit].test(temp_node)||degInG[temp_node]<temp_stand1-k){
                    Cand2->remove(temp_node);
                    remove_C2.push_back(temp_node);
                    for(temp_j=degree[temp_node]-1;temp_j>=0;--temp_j){
                        temp_node2=Graph[temp_node][temp_j];
                        --degInG[temp_node2];
                    }
                }
            }

            for(temp_i=Exc2->vnum-1;temp_i>=0;--temp_i){
                temp_node=Exc2->vlist[temp_i];
                if(!bit_G[povit].test(temp_node)||degInG[temp_node]<temp_stand1-k){
                    Exc2->remove(temp_node);
                    remove_E2.push_back(temp_node);
                }
            }
        }else{
            temp_size1=P1->vnum-k;
            for(temp_i=Cand2->vnum-1;temp_i>=0;--temp_i){
                temp_node=Cand2->vlist[temp_i];
                if(degInP[temp_node]<temp_size1||degInG[temp_node]<temp_stand1-k){
                    Cand2->remove(temp_node);
                    remove_C2.push_back(temp_node);
                    for(temp_j=degree[temp_node]-1;temp_j>=0;--temp_j){
                        temp_node2=Graph[temp_node][temp_j];
                        --degInG[temp_node2];
                    }
                }
            }

            for(temp_i=Exc2->vnum-1;temp_i>=0;--temp_i){
                temp_node=Exc2->vlist[temp_i];
                if(degInP[temp_node]<temp_size1||degInG[temp_node]<temp_stand1-k){
                    Exc2->remove(temp_node);
                    remove_E2.push_back(temp_node);
                }
            }
        }

        
        for(temp_i=degree[povit]-1;temp_i>=0;--temp_i){
            temp_node=Graph[povit][temp_i];
            if(Cand2->contains(temp_node)||P2->contains(temp_node))
                G_index[temp_node]=1;
        }

        temp_vector.clear();
        temp_size1=P1->vnum;
        temp_size2=P2->vnum;
        for(temp_i=P2->vnum-1;temp_i>=0;--temp_i){
            temp_node=P2->vlist[temp_i];
            if(degInP[temp_node]==temp_size1-k&&!bit_G[povit].test(temp_node)){
                temp_vector.push_back(temp_node);
            }
        }

        if(!temp_vector.empty()&&Cand2->vnum+P2->vnum>=temp_stand2){
            for(temp_i=Cand1->vnum-1;temp_i>=0&&Cand1->vnum+P1->vnum>=temp_stand1;--temp_i){
                temp_node=Cand1->vlist[temp_i];
                if(degInP[temp_node]==temp_size2) continue;
                if(degInP[temp_node]<temp_size2-k||degInG[temp_node]<temp_stand2-k){
                    Cand1->remove(temp_node);
                    remove_C1.push_back(temp_node);
                    for(temp_j=degree[temp_node]-1;temp_j>=0;--temp_j){
                        temp_node2=Graph[temp_node][temp_j];
                        --degInG[temp_node2];
                    }
                }else{
                    for(int i:temp_vector){
                        if(!bit_G[temp_node].test(i)){
                            Cand1->remove(temp_node);
                            remove_C1.push_back(temp_node);
                            for(temp_j=degree[temp_node]-1;temp_j>=0;--temp_j){
                                temp_node2=Graph[temp_node][temp_j];
                                --degInG[temp_node2];
                            }
                            break;
                        }
                    }
                    if(Cand1->contains(temp_node)){
                        count=0;
                        for(temp_j=degree[temp_node]-1;temp_j>=0;--temp_j){
                            if(G_index[Graph[temp_node][temp_j]]){
                                count++;
                            }
                        }
                        if(count<temp_stand2-2*k){
                            Cand1->remove(temp_node);
                            remove_C1.push_back(temp_node);
                            for(temp_j=degree[temp_node]-1;temp_j>=0;--temp_j){
                                temp_node2=Graph[temp_node][temp_j];
                                --degInG[temp_node2];
                            }
                        }
                    }
                }
            }

            for(temp_i=Exc1->vnum-1;temp_i>=0;--temp_i){
                temp_node=Exc1->vlist[temp_i];
                if(degInP[temp_node]==temp_size2) continue;
                if(degInP[temp_node]<temp_size2-k||degInG[temp_node]<temp_stand2-k){
                    Exc1->remove(temp_node);
                    remove_E1.push_back(temp_node);
                }else{
                    for(int i:temp_vector){
                        if(!bit_G[temp_node].test(i)){
                            Exc1->remove(temp_node);
                            remove_E1.push_back(temp_node);
                            break;
                        }
                    }
                }
                
            }
        }
        for(temp_i=degree[povit]-1;temp_i>=0;--temp_i){
            G_index[Graph[povit][temp_i]]=0;
        }
       
        if(CanBranch()&&P_R.vnum+Cand_R.vnum>=right&&P_L.vnum+Cand_L.vnum>=left)
            BKmb_Rec();

        for(int i:remove_C1){
            Cand1->add(i);
            for(temp_i=degree[i]-1;temp_i>=0;--temp_i){
                ++degInG[Graph[i][temp_i]];
            }
        }

        for(int i:remove_C2){
            Cand2->add(i);
            for(temp_i=degree[i]-1;temp_i>=0;--temp_i){
                ++degInG[Graph[i][temp_i]];
            }
        }

        P1->remove(povit);
        Cand1->add(povit);
        for(temp_i=degree[povit]-1;temp_i>=0;--temp_i){
            --degInP[Graph[povit][temp_i]];
        }

        for(int i:remove_E1){
            Exc1->add(i);
        }
        for(int i:remove_E2){
            Exc2->add(i);
        }

        

    }

}

bool BKMB::CanBranch(){
    temp_j=P_L.vnum+Cand_L.vnum;
    temp_k=P_R.vnum+Cand_R.vnum;
    for(temp_i=Exc_L.vnum-1;temp_i>=0;--temp_i){
        temp_node=Exc_L.vlist[temp_i];
        if(degInG[temp_node]==temp_k){
            return false;
        }
    }

    for(temp_i=Exc_R.vnum-1;temp_i>=0;--temp_i){
        temp_node=Exc_R.vlist[temp_i];
        if(degInG[temp_node]==temp_j){
            return false;
        }
    }

    return true;
}