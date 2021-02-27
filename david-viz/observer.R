source("load_data.R", local = T)
observe({
  
  temp.mybindeddf = mybindeddf
  updatePickerInput(
    session = session,
    inputId = "SetGenePicker",
    choices = data.frame(
      mybindeddf %>% select(Gene_Name) %>% unique()
    )$Gene_Name
  )
  updatePickerInput(
    session = session,
    inputId = "SetThemePicker",
    choices = data.frame(
      mybindeddf %>% select(Group) %>% unique()
    )$Group
  )
  
})
