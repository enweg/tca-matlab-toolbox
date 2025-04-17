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

Y = SVAR.simulate(T, A0, APlus, 'trendExponents', trendExponents);
model = SVAR(Y, p, 'trendExponents', trendExponents);

model.isFitted()
[A0, APlus] = model.coeffs();
model.fitted();
model.residuals();
model.makeCompanionMatrix();
model.spectralRadius();
model.isStable();
model.aic();
model.bic();
model.sic();
model.hqc();

model.getDependent()
model.getIndependent()
model.getInputData()
model.getVariableNames()
model.nobs() == T - p
model.isStructural()

method = Recursive();
model.fit(method)
model.isFitted()
[A0, APlus] = model.coeffs()
model.fitted()
model.residuals()
model.makeCompanionMatrix()
model.spectralRadius()
model.isStable()
model.aic()
model.bic()
model.sic()
model.hqc()

modelBig = SVAR(Y, p + 10, 'trendExponents', trendExponents);
[modelBest, icTable] = modelBig.fitAndSelect(method, @VAR.aic_)
modelBest.p
icTable

model = SVAR(Y, p, 'trendExponents', trendExponents);
method = Recursive();
model.fit(method);

maxHorizon = 10;
obj = model;
irfObj = model.IRF(10);
irfObj

model = VAR(Y, p, 'trendExponents', trendExponents);
model.fit()
method = Recursive();
irfObj = model.IRF(10, 'identificationMethod', method);
irfObj
method = InternalInstrument(2);
irfObj = model.IRF(10, 'identificationMethod', method);
irfObj
method = InternalInstrument('Y2');
irfObj = model.IRF(10, 'identificationMethod', method);
irfObj
method = InternalInstrument('Y2', 'normalisingHorizon', 1);
irfObj = model.IRF(10, 'identificationMethod', method);
irfObj

