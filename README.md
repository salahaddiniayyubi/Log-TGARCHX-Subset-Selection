# Log-TGARCHX Model Enhanced With Best Subset Selection Procedures
Authors: **Francois Septier & Audrey Poterie & Victor Elvira & Samir Orujov & Farid Rajabov** 

## VS - LTGARCHX: A Flexible Subset Selection Approach for Estimation of Log - TGARCHX Models and Its Application to BTC Markets

Please refer to [the paper](https://www.overleaf.com/read/xkpfrfgtdwhv) for theoretical underpinnings and other details.

## Running The Main Script
* The whole relevant R environment is given [here](Log-GARCH.Rdata). The easiest way to examine the results without running the whole code is to load the latter environment and go through [main R script](./main.R) by running chunks of it by choice and/or inspecting R objects in the environment.   
* Alternatively, if one insists on running the script, then, the preprocessed [data](./data/data_stationary.csv) can be directly utilized to run the [main R script](./main.R) and to get the major findings of our paper. However, running the code on personal computers may take a while. 
* The third alternative is to use parallel computing for fast script runtime. The code is optimized for **parallel computing**, and it was run successfully at [UBS cluster](https://cluster-irisa.univ-ubs.fr/wiki/). The slurm jobs script is [here](./slurm.sh).

## Data Collection and Manipulation Process In More Details

- If one wants to inspect the raw data collection process in more details then [this python code](./data/data_download.ipynb) should be used.
- However, you need to supply your own Binance, Quandl and FRED Api keys to run the latter script successfully.
- After downloading the data using the latter script, the data needs to be pre-processed using [this script](./data/preprocessing.ipynb).

## Few Words About Data

1. Based on previous research, the data sources are classified into 5 groups, namely, Blockchain info, Public Opinion, Risks and Uncertainties, Financials and Macroeconomic development.
2. Data on 39 variables in total were collected in daily frequency. 
3. Moreover, time series of BTC price were collected on 5-min frequency to estimate the realized voatility.
4. The first lag of daily, weekly and monthly realized variances in daily frequencies are also used as conditioning variables as in HAR-RV models.  
5. The sample period runs from December 18, 2017 to June 17, 2022 and includes 1643 observations in total.
6. Data Collection Sources: Yahoo Finance, Quandl, FRED, Binance Exchange.

## Robustness Check

1. Rolling forecasts are asserted to be better while evaluating the nested models. Please refer to the paper for details.
2. We used both rolling and recursive forecasts and it turns out that our methodology works well for both forecast schemes.
3. Please refer to [Rolling window scheme output](./rolling.out) and [Recursive window scheme output](./recursive.out) for details.

Happy Reading!
