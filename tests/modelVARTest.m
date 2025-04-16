addpath("./models")

k = 3;
p = 2;
T = 100;
trendExponents = [0];
B = 0.2 * randn(k, k*p + length(trendExponents));

errors = randn(k, T);
Y = VAR.simulate(errors, B, 'trendExponents', trendExponents);

model = VAR(Y, p, 'trendExponents', trendExponents)

model.isFitted()

% following are expected to throw error because VAR is not estimated
model.coeffs()
model.fitted()
model.residuals()
model.ncoeffs()
model.makeCompanionMatrix()
model.spectralRadius()
model.isStable()
model.aic()
model.bic()
model.hqc()
model.sic()
model.IRF(10);

model.nobs()
model.getDependent()
model.getIndependent()
model.getInputData()

model.isStructural()

model.fit()
% the following should now all work
model.coeffs()
model.coeffs(true)
model.fitted()
model.residuals()
model.ncoeffs()
model.makeCompanionMatrix()
model.spectralRadius()
model.isStable()
model.aic()
model.bic()
model.hqc()
model.sic()

modelLarge = VAR(Y, 10, 'trendExponents', trendExponents)
[modelBest, icTable] = modelLarge.fitAndSelect()
modelBest.p
icTable
[modelBest, icTable] = modelLarge.fitAndSelect(@VAR.hqc_)
modelBest.p
icTable
[modelBest, icTable] = modelLarge.fitAndSelect(@VAR.sic_)
modelBest.p
icTable


model = VAR(Y, p, 'trendExponents', trendExponents)
model.fit()
irfs = model.IRF(10);
irfs
