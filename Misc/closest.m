function idx = closest(Array,Value)

  [~,idx] = min(abs(Array-Value));

end

