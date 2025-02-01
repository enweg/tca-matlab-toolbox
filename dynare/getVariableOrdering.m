function ordering=getVariableOrdering(M_, options_, oo_, vars)
    cellordering = M_.endo_names(oo_.dr.order_var(vars));
    ordering = string(cell2mat(cellordering(1)));
    for i=2:length(cellordering)
        ordering = [ordering; string(cell2mat(cellordering(i)))];
    end
end
