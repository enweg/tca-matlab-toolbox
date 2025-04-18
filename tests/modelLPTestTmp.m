
k = 3; 
p = 2;
T = 1000;
trendExponents = [0];
m = length(trendExponents);

A0 = randn(k, k);
A0 = tril(A0);
S = diag(sign(diag(A0)));
A0 = A0 * S;

B = 0.2 * randn(k, k*p + m);
APlus = A0 * B;

Y =  SVAR.simulate(T, A0, APlus, 'trendExponents', trendExponents);

model = LP(Y, 1, p, 0:2);
model.isFitted()
model.nobs()
model.getDependent()
model.getIndependent()
model.getInputData()
model.getVariableNames()
model.isStructural()

% throws errors because not fitted
model.coeffs()
model.fitted()
model.residuals()

method = Recursive();
model.fit(method);
model.fit();

model.coeffs()

irfObj = model.IRF(2)
irfsLP = irfObj.irfs;

modelLarge = LP(Y, 1, 10, 0:2);
[modelBest, icTable] = modelLarge.fitAndSelect(Recursive(), @VAR.aic_);
modelBest.p
icTable




