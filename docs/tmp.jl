function f()
    document = Documenter.Documents.Document(; sitename="SumOfSquares")
    Documenter.Selectors.runner(Documenter.Builder.SetupBuildDirectory, document)
    Documenter.Selectors.runner(Documenter.Builder.ExpandTemplates, document)
    #Documenter.Selectors.runner(Documenter.Builder.RenderDocument, document)
    document
end
