import sys

import pandas as pd

sys.path.append('models')
import get_data
import torch.utils.data as Data
from utils import train_model,show_info,setup_seed,analizeResult
import Config
from get_data import collate_fn,collate_fn_pair
import utils
import torch
import warnings
warnings.filterwarnings('ignore')
from models import RSSE
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
setup_seed()



class Config(object):
    def __init__(self,dataset):



        self.dataset= dataset # chengdu xian porto
        self.metric='edr' # edr erp dtw
        self.model='TMN'   #T2V,km_test_model,T3S,TMN
        self.task='subtra' # subtra tra

        self.query_min_len = 30 # 10 - 100
        self.query_max_len = 89
        self.data_min_len = 90 # 
        self.data_max_len= 300

        self.train_rate=0.4 # 模型参数
        self.val_rate=0.05

        # path
        self.base_path=r'../'
        self.data_path = self.base_path + 'data' + '/' + self.dataset
        self.result_file = self.base_path+"/"
        self.model_path=self.base_path+"/"
        self.score_path = self.base_path + "/"





        # train_parameter
        self.epoch = 100


        self.num_workers = 0
        self.patience = 5
        self.pair=True



        # data_parameter


        self.device=torch.device('cuda' if torch.cuda.is_available() else 'cpu')

        # 正则化有关参数
        self.x_m=0
        self.x_std=0
        self.y_max = 0
        self.y_std = 0




        # model_parameter
        self.batch_size = 512
        self.embedding_dim = 2
        self.num_hiddens = 64
        self.num_layers = 3
        self.bi=2


def process_data():
    print("调用process_data")
    dataset='xian'
    config = Config.Config(dataset)
    config.metric='edr'
    config.model='RSSE'
    config.pair_p=float(0.03)
    
    score_table, t_table, q_table = get_data.get_data(config)

    train_score_table, val_score_table, test_score_table = get_data.split_data_pair(score_table,config)
    print("train_score_table")
    print(train_score_table.columns)
    print(train_score_table.shape)
    print(train_score_table.head(3))
    print("val_score_table")
    print(val_score_table.columns)
    print(val_score_table.shape)
    print(val_score_table.head(3))
    print("test_score_table")
    print(test_score_table.columns)
    print(test_score_table.shape)
    print(test_score_table.head(3))
        
    train_dataset = get_data.get_Dataset_pair(train_score_table, t_table, q_table,task=config.task)
    val_dataset = get_data.get_Dataset_pair(val_score_table, t_table, q_table,task=config.task)
    test_dataset = get_data.get_Dataset(test_score_table, t_table, q_table,task=config.task)
    return train_dataset,val_dataset,test_dataset


if __name__ == '__main__':
    dataset=sys.argv[1]

    config = Config.Config(dataset)
    config.metric=sys.argv[2]
    config.model=sys.argv[3]
    config.pair_p=float(sys.argv[4])
    config.num_hiddens = int(sys.argv[5])
    # ==data and dateset

    show_info(config)


    if config.model=='RSSE':
        model = RSSE.Siamese_RSSE(
            config.embedding_dim,
            config.num_hiddens,
            config.num_layers
        )
        
        process_data()
        
        
        score_table, t_table, q_table = get_data.get_data(config)


        train_score_table, val_score_table, test_score_table = get_data.split_data_pair(score_table,config)
        print("train_score_table")
        print(train_score_table.columns)
        print(train_score_table.shape)
        print(train_score_table.head(3))
        print("val_score_table")
        print(val_score_table.columns)
        print(val_score_table.shape)
        print(val_score_table.head(3))
        print("test_score_table")
        print(test_score_table.columns)
        print(test_score_table.shape)
        print(test_score_table.head(3))
        

        train_dataset = get_data.get_Dataset_pair(train_score_table, t_table, q_table,task=config.task)
        val_dataset = get_data.get_Dataset_pair(val_score_table, t_table, q_table,task=config.task)
        test_dataset = get_data.get_Dataset(test_score_table, t_table, q_table,task=config.task)
        
        
        
        
        train_loader = Data.DataLoader(
            dataset=train_dataset,
            batch_size=config.batch_size,
            shuffle=True,
            num_workers=4,
            collate_fn=collate_fn_pair,
            drop_last=True,
            pin_memory=True

        )
        val_loader = Data.DataLoader(
            dataset=val_dataset,
            batch_size=config.batch_size,
            shuffle=True,
            num_workers=2,
            collate_fn=collate_fn_pair,
            drop_last=True,
            pin_memory=True

        )
        test_dataloader = Data.DataLoader(
            dataset=test_dataset,
            batch_size=config.batch_size,
            shuffle=False,
            num_workers=2,
            collate_fn=collate_fn,
            drop_last=False,
            pin_memory=True

        )












    show_info(config)

    model = train_model(model, train_loader, val_loader, config)
    eva_data=utils.evaluate(model, test_score_table,test_dataloader,config)

    result_list,columns_list= analizeResult(eva_data)
    result_list.append(config.dataset),columns_list.append('dataset')
    result_list.append(config.metric),columns_list.append('metric')
    result_list.append(config.model), columns_list.append('model')
    result_list=pd.DataFrame(result_list).T
    result_list.columns=columns_list
    result_list.to_csv('/root/copy/subtra/5_1_result/score/'+config.dataset+'_'+config.metric+'_'+config.model+'_'+str(config.num_hiddens)+'_'+str(config.pair_p),index=False)

