ttk_add_paraview_plugin(ttkTensorMonotonousSubdivision
  SOURCES
    ${VTKWRAPPER_DIR}/ttkTensorMonotonousSubdivision/ttkTensorMonotonousSubdivision.cpp
  PLUGIN_XML
    TensorMonotonousSubdivision.xml
  LINK
    tensorMonotonousSubdivision
    )

