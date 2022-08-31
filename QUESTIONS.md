# Questions

## Overall
- [ ] Why do the normalization on X before training?
- [ ] In W decode, why is k-sparsity not enforced (instead lda2/10 is used)?
- [ ] In decode() fncs., why is the error tolerance in spams.lasso() set to lambda1=lda2*Xnorm,
      depending on the norm of X (Xnorm = np.linalg.norm(X)**2 / X.shape[1])?
- [ ] Why is training not evaluated on binarised phi/A ?
- [ ] Observed problem: Covariance before is heightened (e.g. also for segmentation errors)
- [ ] Observed problem: Rare events/rarely expressed proteins are not captured accuratelly



## Small Team

### 10.08.22
- [ ] Should we look at differences of markers, or even differences in pos/neg classes?
- [ ] In simulated X, there are very neg values, does that make sense? Should we put them to 0?
- [ ] What is a good marker to analyse resulting X?
- [ ] Transfomation for plotting?
- [ ] What kind of clipping max is done by default in plotCell()?
