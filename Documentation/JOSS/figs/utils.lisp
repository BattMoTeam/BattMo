;; requires that the convert command from ImageMagick is installed on your system.
(cl-loop for file in (list "batterygeometries" "optimxample")
         do (shell-command (concat "convert -density 600 " file ".pdf " file ".png"))
         )
