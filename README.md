# Designing-an-automated-trading-strategy
Based on daily stock market and sentimental data, an automated trading strategy is proposed with the aim of closing operations with profits.

In total, four different strategies are proposed, ranging from most to least profit generated, which are guided by values such as opening prices, pivot points, volumes of positive and negative news, or a combination of all. Using the functions TestStrategy, TestStrategyComp, and RollingTestStrategy, the effectiveness of each method can be verified. 

It can be observed that the most profitable methods are Sent_and_PP3 and Sent_and_PP4, first one being less risky and therefore generating lower profits, while the second one can be made as aggressive as desired, aiming to generate as much profit as desired with a proportional level of risk.
