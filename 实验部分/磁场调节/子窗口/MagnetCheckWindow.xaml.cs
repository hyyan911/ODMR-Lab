﻿using Controls.Charts;
using ODMR_Lab.Python.LbviewHandler;
using ODMR_Lab.Windows;
using ODMR_Lab.位移台部分;
using ODMR_Lab.基本控件;
using ODMR_Lab.实验部分.扫描基方法;
using ODMR_Lab.实验部分.磁场调节;
using System;
using System.Collections.Generic;
using System.Diagnostics.Contracts;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Media.Media3D;
using System.Windows.Shapes;

namespace ODMR_Lab.磁场调节
{
    /// <summary>
    /// MagnetZWindow.xaml 的交互逻辑
    /// </summary>
    public partial class MagnetCheckWindow : Window
    {
        public CWPointObject CWPoint1 { get; set; } = null;

        public CWPointObject CWPoint2 { get; set; } = null;

        private DisplayPage ParentPage { get; set; } = null;

        private NumricChartData1D LineCW1 = new NumricChartData1D() { Name = "CW1", DataAxisType = ChartDataType.Y };
        private NumricChartData1D LineCW2 = new NumricChartData1D() { Name = "CW2", DataAxisType = ChartDataType.Y };
        private NumricChartData1D Freq1 = new NumricChartData1D() { Name = "Freq1", DataAxisType = ChartDataType.X };
        private NumricChartData1D Freq2 = new NumricChartData1D() { Name = "Freq2", DataAxisType = ChartDataType.X };

        public MagnetCheckWindow(string windowtitle, DisplayPage parent)
        {
            InitializeComponent();
            ParentPage = parent;
            Title = windowtitle;
            WindowTitle.Content = windowtitle;
            CodeHelper.WindowResizeHelper helper = new CodeHelper.WindowResizeHelper();
            helper.RegisterWindow(this, dragHeight: 30);
        }

        private void Close(object sender, RoutedEventArgs e)
        {
            Hide();
        }

        private void Minimize(object sender, RoutedEventArgs e)
        {
            WindowState = WindowState.Minimized;
        }
        #region 刷新显示部分
        /// <summary>
        /// 刷新CW点列表
        /// </summary>
        public void UpdatePointData()
        {
            MeasuredPoints.ClearItems();
            if (CWPoint1 == null)
            {
                MeasuredPoints.AddItem(null, "", "", "", "", "");
            }
            else
            {
                MeasuredPoints.AddItem(null, Math.Round(CWPoint1.CW1, 4).ToString(), Math.Round(CWPoint1.CW2, 4).ToString(), Math.Round(CWPoint1.Bp, 4).ToString(), Math.Round(CWPoint1.Bv, 4).ToString(), Math.Round(CWPoint1.B, 4).ToString());
            }

            if (CWPoint2 == null)
            {
                MeasuredPoints.AddItem(null, "", "", "", "", "");
            }
            else
            {
                MeasuredPoints.AddItem(null, Math.Round(CWPoint2.CW1, 4).ToString(), Math.Round(CWPoint2.CW2, 4).ToString(), Math.Round(CWPoint2.Bp, 4).ToString(), Math.Round(CWPoint2.Bv, 4).ToString(), Math.Round(CWPoint2.B, 4).ToString());
            }
        }

        /// <summary>
        /// 刷新计算
        /// </summary>
        public void UpdateCalculate()
        {
            Dispatcher.Invoke(() =>
            {
                try
                {
                    MagnetScanConfigParams Cp = new MagnetScanConfigParams();
                    Cp.ReadFromPage(new FrameworkElement[] { ParentPage });
                    MagnetScanExpParams Pp = new MagnetScanExpParams();
                    Pp.ReadFromPage(new FrameworkElement[] { ParentPage });
                    MagnetAutoScanHelper.CalculatePossibleLocs(Cp, Pp, out double x1, out double y1, out double z1, out double angle1, out double x2, out double y2, out double z2, out double angle2);
                    PredictParams.ClearItems();
                    PredictParams.AddItem(null, x1, y1, z1, angle1);
                    PredictParams.AddItem(null, x2, y2, z2, angle2);
                }
                catch (Exception ex)
                {
                    MessageWindow.ShowTipWindow("参数格式错误，磁场预测计算未完成", this);
                }
            });
        }

        /// <summary>
        /// 刷新图表显示
        /// </summary>
        public void UpdateChartDisplay()
        {
            if (MeasuredPoints.GetSelectedTag() == null) return;
            CWPointObject point = MeasuredPoints.GetSelectedTag() as CWPointObject;

            LineCW1.Data = point.CW1Values;
            LineCW2.Data = point.CW2Values;
            Freq1.Data = point.CW1Freqs;
            Freq2.Data = point.CW2Freqs;

            CWChart.UpdateChartDataPanel();
            CWChart.UpdateChart(true);
        }

        /// <summary>
        /// 刷新所有显示
        /// </summary>
        public void UpdateDisplay()
        {
            UpdatePointData();
            UpdateChartDisplay();
        }
        #endregion

        #region 移动并测量
        private void MovetoP1(object sender, RoutedEventArgs e)
        {
            double x = 0, y = 0, z = 0, angle = 0;
            try
            {
                x = double.Parse(PredictParams.GetCellValue(0, 0).ToString());
                y = double.Parse(PredictParams.GetCellValue(0, 1).ToString());
                z = double.Parse(PredictParams.GetCellValue(0, 2).ToString());
                angle = double.Parse(PredictParams.GetCellValue(0, 3).ToString());
            }
            catch (Exception ex)
            {
                MessageWindow.ShowTipWindow("移动目标位置非法，无法进行移动", MainWindow.Handle);
                return;
            }

            Thread t = new Thread(() =>
            {
                Dispatcher.Invoke(() =>
                {
                    P1Btn.IsEnabled = false;
                    P2Btn.IsEnabled = false;
                });
                CWPointObject point = MoveAndScan(x, y, z, angle, null);
                if (point != null)
                {
                    CWPoint1 = point;
                }
                Dispatcher.Invoke(() =>
                {
                    UpdateDisplay();
                    P1Btn.IsEnabled = true;
                    P2Btn.IsEnabled = true;
                });
            });
            t.Start();
        }

        private void MovetoP2(object sender, RoutedEventArgs e)
        {
            double x = 0, y = 0, z = 0, angle = 0;
            try
            {
                x = double.Parse(PredictParams.GetCellValue(1, 0).ToString());
                y = double.Parse(PredictParams.GetCellValue(1, 1).ToString());
                z = double.Parse(PredictParams.GetCellValue(1, 2).ToString());
                angle = double.Parse(PredictParams.GetCellValue(1, 3).ToString());
            }
            catch (Exception ex)
            {
                MessageWindow.ShowTipWindow("移动目标位置非法，无法进行移动", MainWindow.Handle);
                return;
            }

            Thread t = new Thread(() =>
            {
                Dispatcher.Invoke(() =>
                {
                    P1Btn.IsEnabled = false;
                    P2Btn.IsEnabled = false;
                });
                CWPointObject point = MoveAndScan(x, y, z, angle, null);
                if (point != null)
                {
                    CWPoint2 = point;
                }
                Dispatcher.Invoke(() =>
                {
                    UpdateDisplay();
                    P1Btn.IsEnabled = true;
                    P2Btn.IsEnabled = true;
                });
            });
            t.Start();
        }

        /// <summary>
        /// 移动到目标位置并进行实验
        /// </summary>
        public CWPointObject MoveAndScan(double x, double y, double z, double a, Action StopMethod)
        {
            try
            {
                MagnetScanConfigParams P = new MagnetScanConfigParams();
                P.ReadFromPage(new FrameworkElement[] { ParentPage });
                //移动位移台
                var moverx = DeviceDispatcher.TryGetMoverDevice(P.XRelate.Value, OperationMode.ReadWrite, PartTypes.Magnnet, true, true);
                if (moverx == null)
                {
                    return null;
                }
                var movery = DeviceDispatcher.TryGetMoverDevice(P.YRelate.Value, OperationMode.ReadWrite, PartTypes.Magnnet, true, true);
                if (movery == null)
                {
                    return null;
                }
                var moverz = DeviceDispatcher.TryGetMoverDevice(P.ZRelate.Value, OperationMode.ReadWrite, PartTypes.Magnnet, true, true);
                if (moverz == null)
                {
                    return null;
                }
                var movera = DeviceDispatcher.TryGetMoverDevice(P.ARelate.Value, OperationMode.ReadWrite, PartTypes.Magnnet, true, true);
                if (movera == null)
                {
                    return null;
                }
                ScanHelper.Move(moverx, StopMethod, P.XRangeLo.Value, P.XRangeHi.Value, x, 10000);
                ScanHelper.Move(movery, StopMethod, P.YRangeLo.Value, P.YRangeHi.Value, y, 10000);
                ScanHelper.Move(moverz, StopMethod, P.ZRangeLo.Value, P.ZRangeHi.Value, z, 10000);
                ScanHelper.Move(movera, StopMethod, -150, 150, a, 10000);
                //测量
                LabviewConverter.AutoTrace(out Exception e);
                MagnetAutoScanHelper.TotalCWPeaks2OrException(out List<double> peaks, out List<double> freqs1, out List<double> contracts1, out List<double> freqs2, out List<double> contracts2);

                CWPointObject point = new CWPointObject(0, Math.Min(peaks[0], peaks[1]), Math.Max(peaks[0], peaks[1]), P.D.Value, freqs1, contracts1, freqs2, contracts2);

                return point;
            }
            catch (Exception ex)
            {
                MessageWindow.ShowTipWindow("移动并测量未完成:" + ex.Message, MainWindow.Handle);
                return null;
            }
        }
        #endregion
    }
}
