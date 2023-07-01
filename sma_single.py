# -*- coding: utf-8 -*-
"""
Created on Sun Nov 13 23:14:34 2022

@author: Vermouth
"""
import logging
import torch
import torch.nn as nn
import numpy as np
from collections import OrderedDict
import pandas as pd
class StructuredAutoencoderNet(nn.Module):
    ## Initialize the network
    def __init__(self, p, encoder_config, decoder_config, dropout_rate = 0):

        super().__init__()
        self.p = p
        self.encoder_config = encoder_config
        self.decoder_config = decoder_config

        ## Save linear layer weights in a list
        self.weights_layer = []

        ## Generate encoder layer 
        index = 0
        self.encoder_layer = []
        for i in range(len(self.encoder_config['dimension']) - 1):
            self.encoder_layer.append(("linear" + str(index), nn.Linear(int(self.encoder_config['dimension'][i]), int(self.encoder_config['dimension'][i + 1]))))
            if i != len(self.encoder_config['dimension']) - 2:
                self.encoder_layer.append(("Sigmoid" + str(index), nn.Sigmoid()))
                self.encoder_layer.append(("dropout" + str(index), nn.Dropout(p = dropout_rate)))
            index += 1
        
        for index, layer in enumerate(self.encoder_layer):
            if layer[0] == "linear":
                self.weights_layer.append(torch.nn.Parameter(layer[1].weight))
                self.encoder_layer[index][1].weight = self.weights_layer[-1]

        ## Generate decoder layer
        index = 0
        self.decoder_layer = []

        for i in range(len(self.decoder_config['dimension']) - 1):
            # decoder_layer.append(("relu" + str(index),nn.ReLU()))
            if i != 0:
                self.decoder_layer.append(("dropout" + str(index), nn.Dropout(p = dropout_rate)))

            self.decoder_layer.append(("linear" + str(index), nn.Linear(int(self.decoder_config['dimension'][i]), int(self.decoder_config['dimension'][i + 1]))))
            if i != len(self.decoder_config['dimension']) - 2:
                self.decoder_layer.append(("Sigmoid" + str(index), nn.Sigmoid()))
            index += 1

        ## encoder_net and decoder_net
        self.encoder_net = nn.Sequential(OrderedDict(
          self.encoder_layer
        ))

        self.decoder_net = nn.Sequential(OrderedDict(
          self.decoder_layer
        ))
    
    # encode and decode function
    def encode(self, X):
        index = 0
        for layer in self.encoder_layer:
            if layer[0] == "linear":
                X = torch.nn.functional.linear(X, self.weights_layer[index])
                index += 1
            else:
                X = layer[1](X)
        # for layer in self.encoder_net:
        #     X = layer(X)
        return X 

    def decode(self, X):
        index = len(self.weights_layer) - 1
        for layer in self.decoder_layer:
            if layer[0] == "linear":
                X = torch.nn.functional.linear(X, self.weights_layer[index].t())
                index -= 1
            else:
                X = layer[1](X)
        # for layer in self.decoder_net:
        #     X = layer(X)
        return X

    # forward network
    def forward(self, X):
        X = self.encode(X)
        X = self.decode(X)
        
        return X 
def mse_loss(input_tensor, target_tensor):
    return ((input_tensor - target_tensor)**2).sum() / input_tensor.size()[0]

def StructuredMaskedAutoencoder(dataset, group, comp_num, N, p, encoder_dimension, decoder_dimension, train_epoch = 1000):
    model = StructuredAutoencoderNet(p, {'dimension' : encoder_dimension}, {'dimension' : decoder_dimension})
    print(model)
    logging.basicConfig(filename="training_18241_512_40_2000_3e-05.log",level=logging.INFO)
    # Training configuration setting
    criterion = torch.nn.MSELoss(reduction='sum')
    optimizer = torch.optim.Adam(model.parameters(), lr = 3e-05)
    
    # Generate tensor list
    tensor_data_list = []
    for i in range(N):
        tensor_data_list.append(torch.tensor(dataset[i].astype(np.float32).T))

    if not isinstance(comp_num, list):
        comp_num = [comp_num]
    
    # Training process
    for epoch in range(int(train_epoch)):
        tensor_pred_list = []
        tensor_pred_list.append(model(tensor_data_list[0]))
        tensor_pred_list.append(model(tensor_data_list[1]))
        loss = 0       
        loss += mse_loss(tensor_data_list[0], tensor_pred_list[0])        
       #loss += criterion(tensor_data_list[0], tensor_pred_list[0])
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        print(f"Epoch {epoch+1}: Training loss: {loss.item()}")
        with torch.no_grad():
            val_loss = mse_loss(tensor_data_list[1], tensor_pred_list[1])

            print(f"Validation loss :{val_loss}")
        logging.info(f"Epoch {epoch+1}: Training Loss: {loss.item()}, Validation Loss :{val_loss}")
    # Embedding
    embedding_list = []
    recovered_list = []
    for i in range(N):
        embedding_list.append(model.encode(tensor_data_list[i]).detach().numpy())
        recovered_list.append(model(tensor_data_list[i]).detach().numpy())
    torch.save(model, "inVitro_sc_autoencoder_18241_512_40_2000_3e-05.pt")
    return embedding_list, recovered_list

