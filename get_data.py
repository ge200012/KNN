# data_loader.py

import pandas as pd
import os

def get_data_score(config):
    base_path = config['data_path']
    print('Loading Data!!!')
    q_idx = '_{}_{}'.format(config['query_min_len'], config['query_max_len'])
    d_idx = '_{}_{}'.format(config['data_min_len'], config['data_max_len'])
    score_table_path = os.path.join(base_path, '{}_{}_'.format(config['query_min_len'], config['data_min_len']) + config['metric'])
    #t_table_path = os.path.join(base_path, 'tra_Frame_td') + d_idx + '_N'
    #q_table_path = os.path.join(base_path, 'tra_Frame_tq') + q_idx + '_N'

    score_table = pd.read_csv(score_table_path)
    score_table.columns = ['que_index', 'tra_idx', 's', 'e', 'score1']

    if config['task'] == 'tra':
        score_table['t'] = score_table['e'] - score_table['s']
        score_table = score_table[score_table['t'] > 2].drop(['t'], axis=1)

    #t_table = pd.read_pickle(t_table_path)
    
    
    #q_table = pd.read_pickle(q_table_path)

    
    #t_table=pd.DataFrame(t_table)
    #t_table_=t_table['tra'].tolist()
    #q_table=pd.DataFrame(q_table)
    #print("score_table")
    #print(score_table.head(3))
    #print("t_table")
    #print(t_table.head(3))
    #print("q_table")
    #print(q_table.head(3))
    
    print('Loading Data  finish!!')
    if config['metric'] == 'edr':
        score_table['score1'] = score_table['score1'] / config['query_max_len']

    return score_table

def get_data_scale_t(config):
    base_path = config['data_path']
    print('Loading Data!!!')
    d_idx = 'scale_t/processed_trajectories.pkl'
    t_table_path = os.path.join(base_path, d_idx)
    
    t_table = pd.read_pickle(t_table_path)
    
    t_table = pd.DataFrame(t_table)
    t_table_ = t_table['tra'].tolist()
    
    print("Loading Data finish!")
    return t_table_
      
def get_data_t(config):
    base_path = config['data_path']
    print('Loading Data!!!')
    #q_idx = '_{}_{}'.format(config['query_min_len'], config['query_max_len'])
    d_idx = '_{}_{}'.format(config['data_min_len'], config['data_max_len'])
    #score_table_path = os.path.join(base_path, '{}_{}_'.format(config['query_min_len'], config['data_min_len']) + config['metric'])
    
    
    
    
    t_table_path = os.path.join(base_path, 'tra_Frame_td') + d_idx
    #t_table_path = os.path.join(base_path, 'tra_Frame_td') + d_idx + '_N'
    
    
    
    
    #q_table_path = os.path.join(base_path, 'tra_Frame_tq') + q_idx + '_N'

    #score_table = pd.read_csv(score_table_path)
    #score_table.columns = ['que_index', 'tra_idx', 's', 'e', 'score1']

    #if config['task'] == 'tra':
    #    score_table['t'] = score_table['e'] - score_table['s']
    #    score_table = score_table[score_table['t'] > 2].drop(['t'], axis=1)
    
    t_table = pd.read_pickle(t_table_path)
    

    # with open("/home/KNN/data/cp_data/trajectory_data",'r') as file:
    #     for i in range(300000):
    #         line = file.readline()
    #         if not line:  # 如果文件行数不足 30
    #             break
    #         print(line.strip())
    # exit(1)
 
    #t_table=pd.read_pickle("/home/KNN/data/cp_data/trajectory_data")
    #q_table = pd.read_pickle(q_table_path)


    t_table=pd.DataFrame(t_table)
    t_table_=t_table['tra'].tolist()
    #q_table=pd.DataFrame(q_table)
    #print("score_table")
    #print(score_table.head(3))
    #print("t_table")
    #print(t_table.head(3))
    #print("q_table")
    #print(q_table.head(3))
    
    print('Loading Data  finish!!')
    #if config['metric'] == 'edr':
    #    score_table['score1'] = score_table['score1'] / config['query_max_len']

    return t_table_

def get_data_q(config):
    base_path = config['data_path']
    print('Loading Data!!!')
    q_idx = '_{}_{}'.format(config['query_min_len'], config['query_max_len'])
    #d_idx = '_{}_{}'.format(config['data_min_len'], config['data_max_len'])
    #score_table_path = os.path.join(base_path, '{}_{}_'.format(config['query_min_len'], config['data_min_len']) + config['metric'])
    #t_table_path = os.path.join(base_path, 'tra_Frame_td') + d_idx + '_N'
    
    
    
    
    q_table_path = os.path.join(base_path, 'tra_Frame_tq') + q_idx 
    #q_table_path = os.path.join(base_path, 'tra_Frame_tq') + q_idx + '_N'
    
    
    
    
    #score_table = pd.read_csv(score_table_path)
    #score_table.columns = ['que_index', 'tra_idx', 's', 'e', 'score1']

    #if config['task'] == 'tra':
    #    score_table['t'] = score_table['e'] - score_table['s']
    #    score_table = score_table[score_table['t'] > 2].drop(['t'], axis=1)

    #t_table = pd.read_pickle(t_table_path)
    
    
    q_table = pd.read_pickle(q_table_path)

    
    q_table=pd.DataFrame(q_table)
    q_table_=q_table['tra'].tolist()
    #q_table=pd.DataFrame(q_table)
    #print("score_table")
    #print(score_table.head(3))
    #print("q_table")
    #print(q_table.head(3))
    #print("q_table")
    #print(q_table.head(3))
    
    print('Loading Data  finish!!')
    #if config['metric'] == 'edr':
    #    score_table['score1'] = score_table['score1'] / config['query_max_len']

    return q_table_