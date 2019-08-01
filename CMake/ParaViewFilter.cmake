# register a new filter to build in the TTK plugin
# also deduce the location of the corresonding vtk.module file
# TODO list of xml file?
macro(ttk_register_pv_filter vtkModule xmlFile)
  if(NOT EXISTS ${VTKWRAPPER_DIR}/${vtkModule}/vtk.module)
    message(FATAL_ERROR
      "Register a paraview filter without the corresponding vtk.module: "
      ${VTKWRAPPER_DIR}/${vtkModule}
      )
  endif()
  list(APPEND TTK_XMLS ${CMAKE_CURRENT_LIST_DIR}/${xmlFile})
  list(APPEND TTK_MODULES ${vtkModule})
  list(APPEND TTK_MODULE_FILES ${VTKWRAPPER_DIR}/${vtkModule}/vtk.module)
endmacro()
