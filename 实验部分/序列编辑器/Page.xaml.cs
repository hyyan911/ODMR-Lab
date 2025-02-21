﻿using CodeHelper;
using Controls;
using Controls.Windows;
using HardWares.仪器列表.板卡.Spincore_PulseBlaster;
using HardWares.温度控制器;
using HardWares.温度控制器.SRS_PTC10;
using HardWares.端口基类;
using HardWares.端口基类部分;
using HardWares.纳米位移台;
using HardWares.纳米位移台.PI;
using ODMR_Lab.Windows;
using ODMR_Lab.基本窗口;
using ODMR_Lab.实验部分.序列编辑器;
using System;
using System.Collections.Generic;
using System.IO;
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
using System.Windows.Navigation;
using System.Windows.Shapes;
using ContextMenu = Controls.ContextMenu;
using Path = System.IO.Path;

namespace ODMR_Lab.序列编辑器
{
    /// <summary>
    /// Page1.xaml 的交互逻辑
    /// </summary>
    public partial class DisplayPage : ExpPageBase
    {
        public override string PageName { get; set; } = "序列编辑器";
        public DisplayPage()
        {
            InitializeComponent();
        }

        public override void InnerInit()
        {
        }

        public override void CloseBehaviour()
        {
        }

        public override void UpdateParam()
        {
        }

        #region 操作
        private void ResizePlot(object sender, RoutedEventArgs e)
        {
            chart.RefreshPlotWithAutoScale();
        }

        private void Snap(object sender, RoutedEventArgs e)
        {
            Clipboard.SetImage(CodeHelper.SnapHelper.GetControlSnap(chart));
            TimeWindow window = new TimeWindow();
            window.Owner = Window.GetWindow(this);
            window.WindowStartupLocation = WindowStartupLocation.CenterOwner;
            window.ShowWindow("截图已复制到剪切板");
        }
        #endregion

        private string GetChannelDescription(SequenceChannel ChannelInd)
        {
            string desc = "Undefined";
            if (MainWindow.Dev_PBPage.PBs.Count != 0)
            {
                desc = MainWindow.Dev_PBPage.PBs[0].FindDescriptionOfChannel((int)ChannelInd);
                if (desc == "")
                    desc = "Undefined";
            }

            return desc;
        }

        SequenceDataAssemble Sequence { get; set; } = new SequenceDataAssemble();
        public void UpdateSequenceData()
        {
            try
            {
                #region 刷新通道列表
                ChannelPanel.ClearItems();
                foreach (var item in Sequence.Channels)
                {

                    var desc = GetChannelDescription(item.ChannelInd);
                    ChannelPanel.AddItem(item, item.ChannelInd, desc, item.IsDisplay);
                }
                #endregion
                SignalPanel.ClearItems();
                SequenceName.Text = Sequence.Name;
                LoopStep.Text = Sequence.LoopCount.ToString();
                RefreshPlot(int.Parse(LoopStep.Text));
            }
            catch { }
        }

        private void ChannelPanel_ItemSelected(int arg1, object arg2)
        {
            SignalPanel.ClearItems();
            SequenceChannelData data = arg2 as SequenceChannelData;
            foreach (var item in data.Peaks)
            {
                SignalPanel.AddItem(item, item.PeakName, item.WaveValue, item.PeakSpan, item.Step);
            }
        }

        private void RefreshPlot(int LoopIndex)
        {
            if (Sequence == null) return;

            List<int> SeqTimes = new List<int>();

            foreach (var item in Sequence.Channels)
            {
                item.ChannelWaveData.Name = Enum.GetName(item.ChannelInd.GetType(), item.ChannelInd);
                item.ChannelWaveData.X.Clear();
                item.ChannelWaveData.Y.Clear();
                if (item.Peaks.Count == 0) continue;
                int timex = 0;
                foreach (var peak in item.Peaks)
                {
                    if (peak.WaveValue == WaveValues.Zero)
                    {
                        item.ChannelWaveData.X.Add(timex);
                        item.ChannelWaveData.Y.Add(0);
                        item.ChannelWaveData.X.Add(timex + peak.PeakSpan + peak.Step * LoopIndex);
                        item.ChannelWaveData.Y.Add(0);
                    }
                    if (peak.WaveValue == WaveValues.One)
                    {
                        item.ChannelWaveData.X.Add(timex);
                        item.ChannelWaveData.Y.Add(1);
                        item.ChannelWaveData.X.Add(timex + peak.PeakSpan + peak.Step * LoopIndex);
                        item.ChannelWaveData.Y.Add(1);
                    }
                    timex += peak.PeakSpan + peak.Step * LoopIndex;
                }
                SeqTimes.Add(timex);
            }
            if (SeqTimes.Count == 0)
            {
                chart.DataList.Clear();
                chart.RefreshPlotWithAutoScale();
                return;
            }
            int maxtime = SeqTimes.Max();
            foreach (var item in Sequence.Channels)
            {
                if (item.ChannelWaveData.GetXCount() == 0) continue;
                if (item.ChannelWaveData.X.Last() <= maxtime)
                {
                    item.ChannelWaveData.X.Add(item.ChannelWaveData.X.Last());
                    item.ChannelWaveData.Y.Add(0);
                    item.ChannelWaveData.X.Add(maxtime);
                    item.ChannelWaveData.Y.Add(0);
                }
            }

            chart.DataList.Clear();
            //重新排列
            double ind = 0;
            foreach (var item in Sequence.Channels)
            {
                if (item.IsDisplay == false) continue;
                item.ChannelWaveData.Y = item.ChannelWaveData.Y.Select(x => x + ind).ToList();
                ind += 1.1;
                chart.DataList.Add(item.ChannelWaveData);
            }
            chart.RefreshPlotWithAutoScale();
        }

        public void UpdatePeakData()
        {
            try
            {
                #region 刷新波形列表
                SignalPanel.ClearItems();
                if (ChannelPanel.GetSelectedTag() == null) return;
                foreach (var seq in (ChannelPanel.GetSelectedTag() as SequenceChannelData).Peaks)
                {
                    SignalPanel.AddItem(seq, seq.PeakName, seq.WaveValue, seq.PeakSpan, seq.Step);
                }
                RefreshPlot(int.Parse(CurrentLoopIndex.Text));
                #endregion
            }
            catch (Exception) { }
        }

        #region 新建波形和通道
        private void NewChannel(object sender, RoutedEventArgs e)
        {
            Sequence.Channels.Add(new SequenceChannelData(SequenceChannel.None));
            UpdateSequenceData();
        }

        private void NewPeak(object sender, RoutedEventArgs e)
        {
            if (ChannelPanel.GetSelectedTag() == null) return;
            var channel = ChannelPanel.GetSelectedTag() as SequenceChannelData;
            channel.Peaks.Add(new SequenceWaveSeg("newseg", 0, 0, WaveValues.Zero, channel));
            UpdatePeakData();
        }
        #endregion

        private void CurrentLoopIndex_KeyDown(object sender, KeyEventArgs e)
        {
            if (e.Key == Key.Enter)
            {
                try
                {
                    RefreshPlot(int.Parse(CurrentLoopIndex.Text));
                }
                catch (Exception)
                {
                }
            }
        }

        private void NewSequence(object sender, RoutedEventArgs e)
        {
            Sequence = new SequenceDataAssemble() { Name = "Newseq" };
            UpdateSequenceData();
            SequencePanel.Visibility = Visibility.Visible;
            TipPanel.Visibility = Visibility.Hidden;
        }

        private void SaveSequence(object sender, RoutedEventArgs e)
        {
            try
            {
                if (SequenceName.Text == "" || SequenceName.Text == null) throw new Exception("序列名不能为空");
                //检查是否已存在同名文件
                DirectoryInfo info = Directory.CreateDirectory(Path.Combine(Environment.CurrentDirectory, "Sequences"));
                var files = info.GetFiles();
                foreach (var item in files)
                {
                    try
                    {
                        var dic = FileObject.ReadDescription(item.FullName);
                        if (dic.Keys.Contains("SequenceAssembleName") && dic["SequenceAssembleName"] == SequenceName.Text)
                        {
                            if (MessageWindow.ShowMessageBox("提示", "已存在同名序列，是否覆盖?", MessageBoxButton.YesNo, owner: Window.GetWindow(this)) == MessageBoxResult.No)
                            {
                                return;
                            }
                            break;
                        }
                    }
                    catch (Exception) { }
                }
                //读取参数
                SequenceDataAssemble a = new SequenceDataAssemble();
                a.Name = SequenceName.Text;
                a.LoopCount = int.Parse(LoopStep.Text);
                foreach (var item in Sequence.Channels)
                {
                    a.Channels.Add(item);
                }
                a.WriteToFile();

                TipPanel.Visibility = Visibility.Visible;
                SequencePanel.Visibility = Visibility.Hidden;

                TimeWindow win = new TimeWindow();
                win.Owner = Window.GetWindow(this);
                win.WindowStartupLocation = WindowStartupLocation.CenterOwner;
                win.ShowWindow("序列已保存");
            }
            catch (Exception ex)
            {
                MessageWindow.ShowTipWindow("保存未完成:" + ex.Message, Window.GetWindow(this));
            }
        }

        private void OpenSequence(object sender, RoutedEventArgs e)
        {
            //保存到目标文件夹
            DirectoryInfo info = Directory.CreateDirectory(Path.Combine(Environment.CurrentDirectory, "Sequences"));
            var files = info.GetFiles();
            List<KeyValuePair<string, string>> values = new List<KeyValuePair<string, string>>();
            foreach (var item in files)
            {
                try
                {
                    var dic = FileObject.ReadDescription(item.FullName);
                    if (dic.Keys.Contains("SequenceAssembleName"))
                    {
                        values.Add(new KeyValuePair<string, string>(dic["SequenceAssembleName"], item.FullName));
                    }
                }
                catch (Exception) { }
            }
            SequenceFileWindow window = new SequenceFileWindow(values);
            window.Owner = Window.GetWindow(this);
            string path = window.ShowDialog();
            if (path == "") return;
            try
            {
                Sequence = SequenceDataAssemble.ReadFromFile(FileObject.ReadFromFile(path));
                UpdateSequenceData();
                SequenceFileName.Content = Sequence.Name;
            }
            catch (Exception) { }
            TipPanel.Visibility = Visibility.Collapsed;
            SequencePanel.Visibility = Visibility.Visible;
        }

        private void WaveFormValueChanged(int arg1, int arg2, object arg3)
        {
            try
            {
                SequenceWaveSeg data = SignalPanel.GetTag(arg1) as SequenceWaveSeg;
                data.PeakName = SignalPanel.GetCellValue(arg1, 0) as string;
                data.WaveValue = (WaveValues)Enum.Parse(typeof(WaveValues), SignalPanel.GetCellValue(arg1, 1) as string);
                data.PeakSpan = int.Parse(SignalPanel.GetCellValue(arg1, 2) as string);
                data.Step = int.Parse(SignalPanel.GetCellValue(arg1, 3) as string);
                RefreshPlot(int.Parse(CurrentLoopIndex.Text));
            }
            catch
            {
            }
        }

        private void ChannelPanel_ItemValueChanged(int arg1, int arg2, object arg3)
        {
            try
            {
                var channel = ChannelPanel.GetTag(arg1) as SequenceChannelData;
                channel.ChannelInd = (SequenceChannel)Enum.Parse(typeof(SequenceChannel), ChannelPanel.GetCellValue(arg1, 0) as string);
                channel.IsDisplay = (bool)ChannelPanel.GetCellValue(arg1, 2);
                ChannelPanel.SetCelValue(arg1, 1, GetChannelDescription(channel.ChannelInd));
                RefreshPlot(int.Parse(CurrentLoopIndex.Text));
            }
            catch (Exception) { }
        }

        private void ChannelPanel_ItemContextMenuSelected(int arg1, int arg2, object arg3)
        {
            //删除
            if (arg1 == 0)
            {
                var seg = arg3 as SequenceWaveSeg;
                seg.ParentChannel.Peaks.Remove(seg);
                UpdateSequenceData();
            }
        }

        private void SignalPanel_ItemContextMenuSelected(int arg1, int arg2, object arg3)
        {
            //删除
            if (arg1 == 0)
            {
                var seg = arg3 as SequenceWaveSeg;
                seg.ParentChannel.Peaks.Remove(seg);
                UpdatePeakData();
            }
            //上方插入
            if (arg1 == 1)
            {
                var seg = arg3 as SequenceWaveSeg;
                seg.ParentChannel.Peaks.Insert(arg1 - 1 < 0 ? 0 : arg1 - 1, new SequenceWaveSeg("newseg", 0, 0, WaveValues.Zero, seg.ParentChannel));
                UpdatePeakData();
            }
            //下方插入
            if (arg1 == 2)
            {
                var seg = arg3 as SequenceWaveSeg;
                seg.ParentChannel.Peaks.Insert(arg1, new SequenceWaveSeg("newseg", 0, 0, WaveValues.Zero, seg.ParentChannel));
                UpdatePeakData();
            }
        }

        private void ChannelPanel_ItemValueChanged_1(int arg1, int arg2, object arg3)
        {

        }
    }
}
