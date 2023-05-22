function Rand = randinrange(Range)


  Rand = rand([1,1]).*range(Range)+min(Range);

end
