# estimacion
Estimacion mediante técnica de MHE
Primero se utilizaba el estimador original para estimar CD a partir de 3 casos: cd cte, cd como función de #M y CD como interpolación de != #M. Luego se le incorporó la estimación de CL. Surge una cuestión con CD porque estima CD vs M y CD vs t y eso no es caractístico del proyectil ya que si varía el ángulo durante la trayectoría también variará el CD. Por lo que se debe estimar CD0 y CDd2. Para esto se modifica las ecuaciones del estimador.